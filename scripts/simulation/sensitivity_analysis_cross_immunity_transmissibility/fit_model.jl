fit_id = isempty(ARGS) ? 241 : parse(Int64, ARGS[1])

using Revise
using DrWatson
using Turing
using DifferentialEquations
using LogExpFunctions
using Random
using CSV
using DataFrames
using Interpolations
using LinearAlgebra
using FillArrays
using immunity_semi_parametric_model

experiment_name = "sensitivity_analysis_cross_immunity_transmissibility"
datadir(args...) = projectdir("data", "simulation", experiment_name, args...)
simulationdir(args...) = projectdir("scripts", "simulation", experiment_name, args...)
results_dir(args...) = projectdir("results", "simulation", experiment_name, args...)
mkpath(results_dir())

model_dict = subset(CSV.read("scripts/simulation/sensitivity_analysis_cross_immunity_transmissibility/model_table.csv", DataFrame; stringtype=String), :fit_id => ByRow(x -> x == fit_id))[1, :] |> x -> Dict(names(x) .=> values(x))

sim_id = model_dict["sim_id"]
max_t = model_dict["max_t"]

prior_takeover_speed = model_dict["prior_takeover_speed"]
data_takeover_speed = model_dict["data_takeover_speed"]

R₀_model = model_dict["R₀_model"]
immunity_model = model_dict["immunity_model"]
CDR_model = model_dict["CDR_model"]

include(simulationdir("shared_constants.txt"))
include(simulationdir("computed_shared_constants.txt"))

include(simulationdir("load_simulated_data.jl"))

include(simulationdir("prior_constants.jl"))
include(projectdir("src", "ode_models", "seirs_log_ode.jl"))
include(projectdir("src", "turing_models", "seirs_super_model_hifi_cdr.jl"))

if data_takeover_speed == "slow"
    logistic_growth_time_offset = logistic_growth_time_offset_slow
elseif data_takeover_speed == "medium"
    logistic_growth_time_offset = logistic_growth_time_offset_medium
elseif data_takeover_speed == "fast"
    logistic_growth_time_offset = logistic_growth_time_offset_fast
end

model_sample = seirs_super_model_hifi_cdr(prob, logistic_growth_time_offset, data_new_cases, data_new_deaths, data_hospitalizations, data_icu, data_new_seq_variant_1, data_new_seq_variant_2, immunity_model, CDR_model, R₀_model, obstimes, param_change_times, seq_obstimes, Tsit5(), 1e-11, 1e-8)

n_samples = 250
n_chains = 4

Random.seed!(1)
prior_samples = sample(model_sample, Prior(), MCMCThreads(), n_samples, n_chains)

mkpath(results_dir("prior_samples_jld2"))
wsave(results_dir("prior_samples_jld2", savename("prior_samples", (@dict fit_id), "jld2")), @dict prior_samples)

Random.seed!(1)
posterior_samples = sample(model_sample, NUTS(-1, 0.8), MCMCThreads(), n_samples, n_chains)

mkpath(results_dir("posterior_samples_jld2"))
mkpath(results_dir("posterior_samples_csv"))

wsave(results_dir("posterior_samples_jld2", savename("posterior_samples", (@dict fit_id), "jld2")), @dict posterior_samples)
CSV.write(results_dir("posterior_samples_csv", savename("posterior_samples", (@dict fit_id), "csv")), DataFrame(posterior_samples))