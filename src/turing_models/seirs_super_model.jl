prob = ODEProblem{true}(seirs_ode_log!,
    zeros(8),
    (0.0, obstimes[end]),
    ones(9))

@model function seirs_super_model(prob, logistic_growth_time_offset, data_new_cases, data_new_deaths, data_hospitalizations, data_icu, data_new_seq_variant_1, data_new_seq_variant_2,
    immunity_model, IHR_model, R₀_model,
    obstimes, param_change_times,
    DEAlgorithm, abstol, reltol)

    l_obstimes = length(obstimes)
    l_param_change_times = length(param_change_times)

    immunity_model ∈ ["constant", "gmrf", "seq-informed"] || throw(ArgumentError("'immunity_model' must be one of 'constant', 'gmrf', or 'seq-informed'"))
    IHR_model ∈ ["constant", "gmrf", "seq-informed"] || throw(ArgumentError("'IHR_model' must be one of 'constant', 'gmrf', or 'seq-informed'"))
    R₀_model ∈ ["constant", "gmrf", "seq-informed"] || throw(ArgumentError("'R₀_model' must be one of 'constant', 'gmrf', or 'seq-informed'"))

    # Priors
    if "seq-informed" ∈ [immunity_model, IHR_model, R₀_model]
        variant_ratio_at_logistic_growth_offset_time_non_centered ~ Normal()
        time_to_saturation_non_centered ~ Normal()
    end

    if immunity_model == "seq-informed"
        prop_variant_2_immunity_offset_non_centered ~ Normal()
    end

    if immunity_model == "constant"
        dur_immune_non_centered ~ Normal()
    elseif immunity_model == "gmrf"
        dur_immune_params_non_centered ~ MvNormal(Zeros(l_param_change_times + 1), I) # +1 for inital value
        σ_dur_immune_non_centered ~ Normal()
    elseif immunity_model == "seq-informed"
        dur_saturated_immune_non_centered ~ Normal()
        dur_immune_shape_non_centered ~ Normal()
        prop_dur_mixed_immune_non_centered ~ Normal()
    end

    if IHR_model == "constant"
        IHR_non_centered ~ Normal()
    elseif IHR_model == "gmrf"
        IHR_params_non_centered ~ MvNormal(Zeros(l_param_change_times + 1), I) # +1 for inital value
        σ_IHR_non_centered ~ Normal()
    elseif IHR_model == "seq-informed"
        IHR_non_centered_variant_1 ~ Normal()
        IHR_non_centered_variant_2 ~ Normal()
    end

    if R₀_model == "constant"
        R₀_non_centered ~ Normal()
    elseif R₀_model == "gmrf"
        R₀_params_non_centered ~ MvNormal(Zeros(l_param_change_times + 1), I) # +1 for inital value
        σ_R₀_non_centered ~ Normal()
    elseif R₀_model == "seq-informed"
        R₀_saturated_non_centered ~ Normal()
        R₀_shape_non_centered ~ Normal()
        prop_R₀_mixed_non_centered ~ Normal()
    end

    dur_latent_non_centered ~ Normal()
    dur_infectious_non_centered ~ Normal()
    dur_hospitalized_non_centered ~ Normal()
    dur_icu_non_centered ~ Normal()
    HICUR_non_centered ~ Normal()
    ICUDR_non_centered ~ Normal()

    ρ_cases_non_centered ~ Normal()
    ρ_deaths_non_centered ~ Normal()
    ρ_seq_non_centered ~ Normal()

    ϕ_hospitalizations_non_centered ~ truncated(Normal(), 0, Inf)
    ϕ_new_cases_non_centered ~ truncated(Normal(), 0, Inf)
    ϕ_new_seq_uncentered ~ truncated(Normal(), 0, Inf)
    ϕ_icu_non_centered ~ truncated(Normal(), 0, Inf)
    ϕ_new_deaths_non_centered ~ truncated(Normal(), 0, Inf)

    # Would rather do these as a DirichletMultinomial
    E_init_prop_non_centered ~ Normal()
    I_init_prop_non_centered ~ Normal()
    R_init_prop_non_centered ~ Normal()

    # Centering tranformations

    dur_latent = exp(dur_latent_non_centered * dur_latent_non_centered_scale + dur_latent_non_centered_loc)
    dur_infectious = exp(dur_infectious_non_centered * dur_infectious_non_centered_scale + dur_infectious_non_centered_loc)
    dur_hospitalized = exp(dur_hospitalized_non_centered * dur_hospitalized_non_centered_scale + dur_hospitalized_non_centered_loc)
    dur_icu = exp(dur_icu_non_centered * dur_icu_non_centered_scale + dur_icu_non_centered_loc)

    ρ_cases = logistic(ρ_cases_non_centered * ρ_cases_non_centered_scale + ρ_cases_non_centered_loc)
    ρ_deaths = logistic(ρ_deaths_non_centered * ρ_deaths_non_centered_scale + ρ_deaths_non_centered_loc)
    ρ_seq = logistic(ρ_seq_non_centered * ρ_seq_non_centered_scale + ρ_seq_non_centered_loc)

    HICUR = logistic(HICUR_non_centered * HICUR_non_centered_scale + HICUR_non_centered_loc)
    ICUDR = logistic(ICUDR_non_centered * ICUDR_non_centered_scale + ICUDR_non_centered_loc)

    ϕ_hospitalizations = ϕ_hospitalizations_non_centered^-2
    ϕ_new_seq = ϕ_new_seq_uncentered^-2
    ϕ_new_cases = ϕ_new_cases_non_centered^-2
    ϕ_icu = ϕ_icu_non_centered^-2
    ϕ_new_deaths = ϕ_new_deaths_non_centered^-2

    E_init_prop = logistic(E_init_prop_non_centered * E_init_prop_non_centered_scale + E_init_prop_non_centered_loc)
    I_init_prop = logistic(I_init_prop_non_centered * I_init_prop_non_centered_scale + I_init_prop_non_centered_loc)
    R_init_prop = logistic(R_init_prop_non_centered * R_init_prop_non_centered_scale + R_init_prop_non_centered_loc)

    if "seq-informed" ∈ [immunity_model, IHR_model, R₀_model]
        # variant2:variant1 ratio
        variant_ratio_at_logistic_growth_offset_time = logistic(variant_ratio_at_logistic_growth_offset_time_non_centered * variant_ratio_at_logistic_growth_offset_time_non_centered_scale + variant_ratio_at_logistic_growth_offset_time_non_centered_loc)
        time_to_saturation = exp(time_to_saturation_non_centered * time_to_saturation_non_centered_scale + time_to_saturation_non_centered_loc)

        logistic_growth_intercept = log(variant_ratio_at_logistic_growth_offset_time)
        logistic_growth_slope = (logit(0.99) - logit(0.01)) / time_to_saturation
        prop_variant_2_fn(x) = logistic.(logistic_growth_intercept .+ logistic_growth_slope .* (x .- logistic_growth_time_offset))
        prop_variant_2 = prop_variant_2_fn(vcat(0, obstimes))
    end

    if immunity_model == "constant"
        dur_immune_init = exp(dur_immune_non_centered * dur_immune_non_centered_scale + dur_immune_non_centered_loc)
        dur_immune_t = fill(dur_immune_init, l_param_change_times + 1)
    elseif immunity_model == "gmrf"
        σ_dur_immune = exp(σ_dur_immune_non_centered * σ_dur_immune_non_centered_scale + σ_dur_immune_non_centered_loc)
        dur_immune_init_non_centered = dur_immune_params_non_centered[1]
        log_dur_immune_steps_non_centered = dur_immune_params_non_centered[2:end]
        dur_immune_init = exp(dur_immune_init_non_centered * dur_immune_init_non_centered_scale + dur_immune_init_non_centered_loc)
        dur_immune_t_no_init = exp.(log(dur_immune_init) .+ cumsum(log_dur_immune_steps_non_centered) * σ_dur_immune)
        dur_immune_t = vcat(dur_immune_init, dur_immune_t_no_init)
    elseif immunity_model == "seq-informed"
        prop_variant_2_immunity_offset = exp(prop_variant_2_immunity_offset_non_centered * prop_variant_2_immunity_offset_non_centered_scale + prop_variant_2_immunity_offset_non_centered_loc)
        prop_variant_2_with_immunity_offset = prop_variant_2_fn(vcat(0, obstimes) .- prop_variant_2_immunity_offset)
        dur_saturated_immune = exp(dur_saturated_immune_non_centered * dur_saturated_immune_non_centered_scale + dur_saturated_immune_non_centered_loc)
        prop_dur_mixed_immune = logistic(prop_dur_mixed_immune_non_centered * prop_dur_mixed_immune_non_centered_scale + prop_dur_mixed_immune_non_centered_loc)
        dur_mixed_immune = dur_saturated_immune * prop_dur_mixed_immune
        dur_immune_shape = exp(dur_immune_shape_non_centered * dur_immune_shape_non_centered_scale + dur_immune_shape_non_centered_loc)

        dur_immune_α₀ = log(dur_saturated_immune)
        dur_immune_α₁ = (log(dur_mixed_immune) - log(dur_saturated_immune)) * 4^dur_immune_shape

        dur_immune_t = exp.(dur_immune_α₀ .+ dur_immune_α₁ .* (prop_variant_2_with_immunity_offset .* (1 .- prop_variant_2_with_immunity_offset)) .^ dur_immune_shape)
    end

    if IHR_model == "constant"
        IHR_init = exp(IHR_non_centered * IHR_non_centered_scale + IHR_non_centered_loc)
        IHR_t = fill(IHR_init, l_param_change_times + 1)
    elseif IHR_model == "gmrf"
        σ_IHR = exp(σ_IHR_non_centered * σ_IHR_non_centered_scale + σ_IHR_non_centered_loc)
        IHR_init_non_centered = IHR_params_non_centered[1]
        logit_IHR_steps_non_centered = IHR_params_non_centered[2:end]
        IHR_init = logistic(IHR_init_non_centered * IHR_init_non_centered_scale + IHR_init_non_centered_loc)
        IHR_t_no_init = logistic.(logit(IHR_init) .+ cumsum(logit_IHR_steps_non_centered) * σ_IHR)
        IHR_t = vcat(IHR_init, IHR_t_no_init)
    elseif IHR_model == "seq-informed"
        prop_variant_2_with_ihr_offset = prop_variant_2_fn(vcat(0, obstimes) .- dur_infectious)
        IHR_variant_1 = logistic(IHR_non_centered_variant_1 * IHR_non_centered_variant_1_scale + IHR_non_centered_variant_1_loc)
        IHR_variant_2 = logistic(IHR_non_centered_variant_2 * IHR_non_centered_variant_2_scale + IHR_non_centered_variant_2_loc)
        IHR_t = (1 .- prop_variant_2_with_ihr_offset) .* IHR_variant_1 .+ prop_variant_2_with_ihr_offset .* IHR_variant_2
    end

    if R₀_model == "constant"
        R₀_init = exp(R₀_non_centered * R₀_non_centered_scale + R₀_non_centered_loc)
        R₀_t = fill(R₀_init, l_param_change_times + 1)

    elseif R₀_model == "gmrf"
        σ_R₀ = exp(σ_R₀_non_centered * σ_R₀_non_centered_scale + σ_R₀_non_centered_loc)
        R₀_init_non_centered = R₀_params_non_centered[1]
        log_R₀_steps_non_centered = R₀_params_non_centered[2:end]
        R₀_init = exp(R₀_init_non_centered * R₀_init_non_centered_scale + R₀_init_non_centered_loc)

        R₀_t_no_init = exp.(log(R₀_init) .+ cumsum(log_R₀_steps_non_centered) * σ_R₀)
        R₀_t = vcat(R₀_init, R₀_t_no_init)

    elseif R₀_model == "seq-informed"
        R₀_saturated = exp(R₀_saturated_non_centered * R₀_saturated_non_centered_scale + R₀_saturated_non_centered_loc)
        prop_R₀_mixed = logistic(prop_R₀_mixed_non_centered * prop_R₀_mixed_non_centered_scale + prop_R₀_mixed_non_centered_loc)
        R₀_mixed = R₀_saturated / prop_R₀_mixed
        R₀_shape = exp(R₀_shape_non_centered * R₀_shape_non_centered_scale + R₀_shape_non_centered_loc)

        R₀_α₀ = log(R₀_saturated)
        R₀_α₁ = (log(R₀_mixed) - log(R₀_saturated)) * 4^R₀_shape

        R₀_t = exp.(R₀_α₀ .+ R₀_α₁ .* (prop_variant_2 .* (1 .- prop_variant_2)) .^ R₀_shape)
    end

    # Natural scale transformation

    γ = 1 / dur_latent
    ν = 1 / dur_infectious
    η = 1 / dur_hospitalized
    ω = 1 / dur_icu

    κ_t = 1 ./ dur_immune_t
    κ_init = κ_t[1]
    κ_t_no_init = κ_t[2:end]

    IHR_init = IHR_t[1]
    IHR_t_no_init = IHR_t[2:end]

    β_t = R₀_t .* ν
    β_init = β_t[1]
    β_t_no_init = β_t[2:end]

    # Initial Condititons
    E_init = E_init_prop * popsize
    I_init = I_init_prop * popsize
    R_init = R_init_prop * popsize
    # H_init is loaded as a constant
    # ICU_init is loaded as a constant
    # D_init is loaded as a constant
    # C_init is loaded as a constant
    # During optimization, E_init, I_init, R_init can sometimes get huge, such that S_init is negative
    # Should reconsider constructing this as a stick-breaking or dirichilet or some other strictly bounded prior
    S_init = max(1, remaining_population_init - (E_init + I_init + R_init))

    u0_reg_scale = [S_init, E_init, I_init, H_init, R_init, C_init, ICU_init, D_init]
    u0 = log.(u0_reg_scale)
    p = [β_init, γ, ν, η, IHR_init, κ_init, HICUR, ω, ICUDR]

    function param_affect_κ!(integrator)
        ind_t = searchsortedfirst(param_change_times, integrator.t) # Find the index of param_change_times that contains the current timestep
        integrator.p[1] = β_t_no_init[ind_t] # Replace β with a new value from β_t_no_init
        integrator.p[5] = IHR_t_no_init[ind_t] # Replace IHR with a new value from IHR_t_no_init
        integrator.p[6] = κ_t_no_init[ind_t] # Replace κ with a new value from κ_t_no_init
    end
    param_callback = PresetTimeCallback(param_change_times, param_affect_κ!, save_positions=(false, false))

    sol = solve(prob, DEAlgorithm; callback=param_callback, saveat=obstimes, save_start=true, verbose=false, abstol=abstol, reltol=reltol, u0=u0, p=p, tspan=(0.0, obstimes[end]))

    if sol.retcode != :Success
        Turing.@addlogprob! -Inf
        return
    end

    # Likelihood calculations
    sol_reg_scale_array = exp.(Array(sol))

    sol_hospitalizations = sol_reg_scale_array[4, 2:end]
    sol_new_cases = sol_reg_scale_array[6, 2:end] - sol_reg_scale_array[6, 1:(end-1)]
    sol_icu = sol_reg_scale_array[7, 2:end]
    sol_new_deaths = sol_reg_scale_array[8, 2:end] - sol_reg_scale_array[8, 1:(end-1)]

    hospitalizations_mean = sol_hospitalizations
    new_cases_mean = sol_new_cases .* ρ_cases
    new_seq_mean = new_cases_mean .* ρ_seq
    icu_mean = sol_icu
    new_deaths_mean = sol_new_deaths .* ρ_deaths

    for i in 1:l_obstimes
        data_hospitalizations[i] ~ NegativeBinomial2(hospitalizations_mean[i], ϕ_hospitalizations)
        data_new_cases[i] ~ NegativeBinomial2(new_cases_mean[i], ϕ_new_cases)
        data_icu[i] ~ NegativeBinomial2(icu_mean[i], ϕ_icu)
        data_new_deaths[i] ~ NegativeBinomial2(new_deaths_mean[i], ϕ_new_deaths)
    end

    if "seq-informed" ∈ [immunity_model, IHR_model, R₀_model]
        new_seq_variant_1_mean = new_seq_mean .* (1 .- prop_variant_2[2:end])
        new_seq_variant_2_mean = new_seq_mean .* prop_variant_2[2:end]

        for i in 1:l_obstimes
            data_new_seq_variant_1[i] ~ NegativeBinomial2(new_seq_variant_1_mean[i], ϕ_new_seq)
            data_new_seq_variant_2[i] ~ NegativeBinomial2(new_seq_variant_2_mean[i], ϕ_new_seq)
        end
    end

    # Generated quantities
    S_compartment = sol_reg_scale_array[1, :]
    E_compartment = sol_reg_scale_array[2, :]
    I_compartment = sol_reg_scale_array[3, :]
    H_compartment = sol_reg_scale_array[4, :]
    R_compartment = sol_reg_scale_array[5, :]
    C_compartment = sol_reg_scale_array[6, :]
    ICU_compartment = sol_reg_scale_array[7, :]
    D_compartment = sol_reg_scale_array[8, :]

    dur_latent_days = dur_latent * time_interval_in_days
    dur_infectious_days = dur_infectious * time_interval_in_days
    dur_hospitalized_days = dur_hospitalized * time_interval_in_days
    dur_immune_days_t = dur_immune_t .* time_interval_in_days
    dur_icu_days = dur_icu * time_interval_in_days

    gq = (
        R₀_t=R₀_t,
        dur_immune_t=dur_immune_t,
        IHR_t=IHR_t,
        ρ_cases=ρ_cases,
        ρ_deaths=ρ_deaths,
        ρ_seq=ρ_seq,
        HICUR=HICUR,
        ICUDR=ICUDR,
        ϕ_hospitalizations=ϕ_hospitalizations,
        ϕ_new_seq=ϕ_new_seq,
        ϕ_new_cases=ϕ_new_cases,
        ϕ_icu=ϕ_icu,
        ϕ_new_deaths=ϕ_new_deaths,
        hospitalizations_mean=hospitalizations_mean,
        new_cases_mean=new_cases_mean,
        new_seq_mean=new_seq_mean,
        icu_mean=icu_mean,
        new_deaths_mean=new_deaths_mean,
        S_compartment=S_compartment,
        E_compartment=E_compartment,
        I_compartment=I_compartment,
        H_compartment=H_compartment,
        R_compartment=R_compartment,
        C_compartment=C_compartment,
        ICU_compartment=ICU_compartment,
        D_compartment=D_compartment,
        dur_latent_days=dur_latent_days,
        dur_infectious_days=dur_infectious_days,
        dur_hospitalized_days=dur_hospitalized_days,
        dur_immune_days_t=dur_immune_days_t,
        dur_icu_days=dur_icu_days
    )


    if "seq-informed" ∈ [immunity_model, IHR_model, R₀_model]
        gq = merge(gq, (
            variant_ratio_at_logistic_growth_offset_time=variant_ratio_at_logistic_growth_offset_time,
            time_to_saturation=time_to_saturation,
            logistic_growth_intercept=logistic_growth_intercept,
            logistic_growth_slope=logistic_growth_slope,
            new_seq_variant_1_mean=new_seq_variant_1_mean,
            new_seq_variant_2_mean=new_seq_variant_2_mean,
            prop_variant_2=prop_variant_2
        ))
    end

    if immunity_model == "gmrf"
        gq = merge(gq, (σ_dur_immune=σ_dur_immune,))
    elseif immunity_model == "seq-informed"
        gq = merge(gq, (
            prop_variant_2_immunity_offset=prop_variant_2_immunity_offset,
            dur_saturated_immune=dur_saturated_immune,
            prop_dur_mixed_immune=prop_dur_mixed_immune,
            dur_mixed_immune=dur_mixed_immune,
            dur_immune_shape=dur_immune_shape,
            dur_immune_α₀=dur_immune_α₀,
            dur_immune_α₁=dur_immune_α₁
        ))
    end

    if IHR_model == "gmrf"
        gq = merge(gq, (σ_IHR=σ_IHR,))
    elseif IHR_model == "seq-informed"
        gq = merge(gq, (
            IHR_variant_1=IHR_variant_1,
            IHR_variant_2=IHR_variant_2
        ))

    end

    if R₀_model == "gmrf"
        gq = merge(gq, (σ_R₀=σ_R₀,))
    elseif R₀_model == "seq-informed"
        gq = merge(gq, (
            R₀_saturated=R₀_saturated,
            prop_R₀_mixed=prop_R₀_mixed,
            R₀_mixed=R₀_mixed,
            R₀_shape=R₀_shape,
            R₀_α₀=R₀_α₀,
            R₀_α₁=R₀_α₁
        ))
    end

    return gq
end