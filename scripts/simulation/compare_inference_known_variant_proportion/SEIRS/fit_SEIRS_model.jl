max_t = isempty(ARGS) ? 28 : parse(Int64, ARGS[1])
using Revise
using DrWatson
using Turing
using DifferentialEquations
using LogExpFunctions
using Random
using CSV
using DataFrames
using Interpolations
using immunity_semi_parametric_model

datadir(args...) = projectdir("data", "simulation", "compare_inference_known_variant_proportion", args...)
simulationdir(args...) = projectdir("scripts", "simulation", "compare_inference_known_variant_proportion", args...)
results_dir(args...) = projectdir("results", "simulation", "compare_inference_known_variant_proportion", "SEIRS", args...)
mkpath(results_dir())

sim_id = 1
include(simulationdir("load_simulated_data.jl"))

include(simulationdir("SEIRS", "prior_constants_SEIRS_model.jl"))
include(projectdir("src", "seirs_log_ode.jl"))
include(projectdir("src", "seirs.jl"))

model_sample = seirs(prob, data_new_cases, data_new_deaths, data_hospitalizations, data_icu, obstimes, Tsit5(), 1e-11, 1e-8)

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