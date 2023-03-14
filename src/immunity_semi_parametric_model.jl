module immunity_semi_parametric_model
using DynamicPPL
using MCMCChains
using Distributions
using ForwardDiff
function NegativeBinomial2(μ, ϕ)
    p = 1 / (1 + μ / ϕ)
    r = ϕ
    # r = clamp(r, nextfloat(zero(r)), prevfloat(typemax(r)))
    # p = clamp(p, nextfloat(zero(p)), one(p))
    # println("μ: ", ForwardDiff.value(μ))
    # println("ϕ: ", ForwardDiff.value(ϕ))
    # println("p: ", ForwardDiff.value(p))
    Distributions.NegativeBinomial(r, p)
end
export NegativeBinomial2


# utilities for working with Turing model parameter names using only the DynamicPPL API

"""
    flattened_varnames_list(model::DynamicPPL.Model) -> Vector{Symbol}

Get a vector of varnames as `Symbol`s with one-to-one correspondence to the
flattened parameter vector.

```julia
julia> @model function demo()
           s ~ Dirac(1)
           x = Matrix{Float64}(undef, 2, 4)
           x[1, 1] ~ Dirac(2)
           x[2, 1] ~ Dirac(3)
           x[3] ~ Dirac(4)
           y ~ Dirac(5)
           x[4] ~ Dirac(6)
           x[:, 3] ~ arraydist([Dirac(7), Dirac(8)])
           x[[2, 1], 4] ~ arraydist([Dirac(9), Dirac(10)])
           return s, x, y
       end
demo (generic function with 2 methods)

julia> flattened_varnames_list(demo())
10-element Vector{Symbol}:
 :s
 Symbol("x[1,1]")
 Symbol("x[2,1]")
 Symbol("x[3]")
 Symbol("x[4]")
 Symbol("x[:,3][1]")
 Symbol("x[:,3][2]")
 Symbol("x[[2, 1],4][1]")
 Symbol("x[[2, 1],4][2]")
 :y
```
"""
function flattened_varnames_list(model::DynamicPPL.Model)
    varnames_ranges = varnames_to_ranges(model)
    nsyms = maximum(maximum, values(varnames_ranges))
    syms = Vector{Symbol}(undef, nsyms)
    for (var_name, range) in varnames_to_ranges(model)
        sym = Symbol(var_name)
        if length(range) == 1
            syms[range[begin]] = sym
            continue
        end
        for i in eachindex(range)
            syms[range[i]] = Symbol("$sym[$i]")
        end
    end
    return syms
end
export flattened_varnames_list

# code snippet shared by @torfjelde
"""
    varnames_to_ranges(model::DynamicPPL.Model)
    varnames_to_ranges(model::DynamicPPL.VarInfo)
    varnames_to_ranges(model::DynamicPPL.Metadata)

Get `Dict` mapping variable names in model to their ranges in a corresponding parameter vector.

# Examples

```julia
julia> @model function demo()
           s ~ Dirac(1)
           x = Matrix{Float64}(undef, 2, 4)
           x[1, 1] ~ Dirac(2)
           x[2, 1] ~ Dirac(3)
           x[3] ~ Dirac(4)
           y ~ Dirac(5)
           x[4] ~ Dirac(6)
           x[:, 3] ~ arraydist([Dirac(7), Dirac(8)])
           x[[2, 1], 4] ~ arraydist([Dirac(9), Dirac(10)])
           return s, x, y
       end
demo (generic function with 2 methods)

julia> demo()()
(1, Any[2.0 4.0 7 10; 3.0 6.0 8 9], 5)

julia> varnames_to_ranges(demo())
Dict{AbstractPPL.VarName, UnitRange{Int64}} with 8 entries:
  s           => 1:1
  x[4]        => 5:5
  x[:,3]      => 6:7
  x[1,1]      => 2:2
  x[2,1]      => 3:3
  x[[2, 1],4] => 8:9
  x[3]        => 4:4
  y           => 10:10
```
"""
function varnames_to_ranges end

varnames_to_ranges(model::DynamicPPL.Model) = varnames_to_ranges(DynamicPPL.VarInfo(model))
varnames_to_ranges(varinfo::DynamicPPL.UntypedVarInfo) = varnames_to_ranges(varinfo.metadata)
function varnames_to_ranges(varinfo::DynamicPPL.TypedVarInfo)
    offset = 0
    dicts = map(varinfo.metadata) do md
        vns2ranges = varnames_to_ranges(md)
        vals = collect(values(vns2ranges))
        vals_offset = map(r -> offset .+ r, vals)
        offset += reduce((curr, r) -> max(curr, r[end]), vals; init=0)
        Dict(zip(keys(vns2ranges), vals_offset))
    end

    return reduce(merge, dicts)
end
function varnames_to_ranges(metadata::DynamicPPL.Metadata)
    idcs = map(Base.Fix1(getindex, metadata.idcs), metadata.vns)
    ranges = metadata.ranges[idcs]
    return Dict(zip(metadata.vns, ranges))
end
export varnames_to_ranges

generate_names(val) = generate_names("", val)
generate_names(vn_str::String, val::Real) = [vn_str;]
function generate_names(vn_str::String, val::NamedTuple)
    return map(keys(val)) do k
        generate_names("$(vn_str)$(k)", val[k])
    end
end
function generate_names(vn_str::String, val::AbstractArray{<:Real})
    results = String[]
    for idx in CartesianIndices(val)
        s = join(idx.I, ",")
        push!(results, "$vn_str[$s]")
    end
    return results
end

function generate_names(vn_str::String, val::AbstractArray{<:AbstractArray})
    results = String[]
    for idx in CartesianIndices(val)
        s1 = join(idx.I, ",")
        inner_results = map(f("", val[idx])) do s2
            "$vn_str[$s1]$s2"
        end
        append!(results, inner_results)
    end
    return results
end

flatten(val::Real) = [val;]
function flatten(val::AbstractArray{<:Real})
    return mapreduce(vcat, CartesianIndices(val)) do i
        val[i]
    end
end
function flatten(val::AbstractArray{<:AbstractArray})
    return mapreduce(vcat, CartesianIndices(val)) do i
        flatten(val[i])
    end
end

function vectup2chainargs(ts::AbstractVector{<:NamedTuple})
    ks = keys(first(ts))
    vns = mapreduce(vcat, ks) do k
        generate_names(string(k), first(ts)[k])
    end
    vals = map(eachindex(ts)) do i
        mapreduce(vcat, ks) do k
            flatten(ts[i][k])
        end
    end
    arr_tmp = reduce(hcat, vals)'
    arr = reshape(arr_tmp, (size(arr_tmp)..., 1)) # treat as 1 chain
    return Array(arr), vns
end

function vectup2chainargs(ts::AbstractMatrix{<:NamedTuple})
    num_samples, num_chains = size(ts)
    res = map(1:num_chains) do chain_idx
        vectup2chainargs(ts[:, chain_idx])
    end

    vals = getindex.(res, 1)
    vns = getindex.(res, 2)

    # Verify that the variable names are indeed the same
    vns_union = reduce(union, vns)
    @assert all(isempty.(setdiff.(vns, Ref(vns_union)))) "variable names differ between chains"

    arr = cat(vals...; dims = 3)

    return arr, first(vns)
end

function MCMCChains.Chains(ts::AbstractArray{<:NamedTuple})
    return MCMCChains.Chains(vectup2chainargs(ts)...)
end

end
