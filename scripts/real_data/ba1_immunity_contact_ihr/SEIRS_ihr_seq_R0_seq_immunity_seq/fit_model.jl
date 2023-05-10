max_t = isempty(ARGS) ? 30 : parse(Int64, ARGS[1])
# need to rewrite to parse county id as well
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
R₀_model = "seq-informed"

datadir(args...) = projectdir("data", "real_data", "ba1", args...)
simulationdir(args...) = projectdir("scripts", "real_data", "ba1_immunity_contact_ihr", args...)
results_dir(args...) = projectdir("results", "real_data", "ba1_immunity_contact_ihr", "SEIRS_ihr_seq_R0_seq_immunity_seq", args...)
mkpath(results_dir())

include(datadir("shared_constants.txt"))
include(simulationdir("load_real_data.jl"))
include(simulationdir("SEIRS_ihr_seq_R0_seq_immunity_seq", "prior_constants.jl"))
include(projectdir("src", "ode_models", "seirs_log_ode.jl"))
include(projectdir("src", "turing_models", "seirs_super_model_hifi.jl"))

model_sample = seirs_super_model_hifi(prob, logistic_growth_time_offset, data_new_cases, data_new_deaths, data_hospitalizations, data_icu, data_new_seq_variant_1, data_new_seq_variant_2, immunity_model, IHR_model, R₀_model, obstimes, param_change_times, seq_obstimes, Tsit5(), 1e-11, 1e-8)

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