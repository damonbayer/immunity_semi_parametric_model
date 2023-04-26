max_t = isempty(ARGS) ? 28 : parse(Int64, ARGS[1])
# max out at 37
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

immunity_model = "seq-informed"
IHR_model = "seq-informed"
R₀_model = "constant"

logistic_growth_time_offset = 23

datadir(args...) = projectdir("data", "simulation", "compare_inference_unknown_variant_proportion_ihr", args...)
simulationdir(args...) = projectdir("scripts", "simulation", "compare_inference_unknown_variant_proportion_ihr", args...)
results_dir(args...) = projectdir("results", "simulation", "compare_inference_unknown_variant_proportion_ihr", "SEIRS_ihr_seq_R0_constant_immunity_seq_flex", args...)
mkpath(results_dir())

sim_id = 1
include(simulationdir("load_simulated_data.jl"))
include(simulationdir("SEIRS_ihr_seq_R0_constant_immunity_seq_flex", "prior_constants.jl"))
include(projectdir("src", "ode_models", "seirs_log_ode.jl"))
include(projectdir("src", "turing_models", "seirs_super_model_flex.jl"))

model_sample = seirs_super_model_flex(prob, logistic_growth_time_offset, data_new_cases, data_new_deaths, data_hospitalizations, data_icu, data_new_seq_variant_1, data_new_seq_variant_2, immunity_model, IHR_model, R₀_model, obstimes, param_change_times, Tsit5(), 1e-11, 1e-8)

n_samples = 250
n_chains = 4

Random.seed!(1)
prior_samples = sample(model_sample, Prior(), MCMCThreads(), n_samples, n_chains)

mkpath(results_dir("prior_samples_jld2"))
wsave(results_dir("prior_samples_jld2", savename("prior_samples", (@dict max_t), "jld2")), @dict prior_samples)

Random.seed!(1)
posterior_samples = sample(model_sample, NUTS(-1, 0.8), MCMCThreads(), n_samples, n_chains)

mkpath(results_dir("posterior_samples_jld2"))
mkpath(results_dir("posterior_samples_csv"))

wsave(results_dir("posterior_samples_jld2", savename("posterior_samples", (@dict max_t), "jld2")), @dict posterior_samples)
CSV.write(results_dir("posterior_samples_csv", savename("posterior_samples", (@dict max_t), "csv")), DataFrame(posterior_samples))