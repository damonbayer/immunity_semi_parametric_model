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

n_forecast_weeks = 4

datadir(args...) = projectdir("data", "simulation", "compare_inference_known_variant_proportion", args...)
simulationdir(args...) = projectdir("scripts", "simulation", "compare_inference_known_variant_proportion", args...)
results_dir(args...) = projectdir("results", "simulation", "compare_inference_known_variant_proportion", "SEIRS", args...)
mkpath(results_dir())

sim_id = 1
include(simulationdir("load_simulated_data.jl"))

include(simulationdir("SEIRS", "prior_constants_SEIRS_model.jl"))
include(projectdir("src", "ode_models", "seirs_log_ode.jl"))
include(projectdir("src", "turing_models", "seirs.jl"))

model_sample = seirs(prob, data_new_cases, data_new_deaths, data_hospitalizations, data_icu, obstimes, Tsit5(), 1e-11, 1e-8)
model_optimization = seirs(prob, data_new_cases, data_new_deaths, data_hospitalizations, data_icu, obstimes, Tsit5(), 1e-13, 1e-10)
model_generated_quantities = seirs(prob, data_new_cases_forecast, data_new_deaths_forecast, data_hospitalizations_forecast, data_icu_forecast, obstimes_forecast, Tsit5(), 1e-13, 1e-10)
model_predict = seirs(prob, missing_data_new_cases_forecast, missing_data_new_deaths_forecast, missing_data_hospitalizations_forecast, missing_data_icu_forecast, obstimes_forecast, Tsit5(), 1e-13, 1e-10)

prior_samples = load(results_dir("prior_samples_jld2", savename("prior_samples", (@dict max_t), "jld2")))["prior_samples"]
posterior_samples = load(results_dir("posterior_samples_jld2", savename("posterior_samples", (@dict max_t), "jld2")))["posterior_samples"]

prior_generated_quantities = Chains(generated_quantities(model_generated_quantities, prior_samples));
mkpath(results_dir("prior_generated_quantities_csv"))
CSV.write(results_dir("prior_generated_quantities_csv", savename("prior_generated_quantities", (@dict max_t), "csv")), DataFrame(prior_generated_quantities))

Random.seed!(1)
prior_predictive = predict(model_predict, prior_samples);
mkpath(results_dir("prior_predictive_csv"))
CSV.write(results_dir("prior_predictive_csv", savename("prior_predictive", (@dict max_t), "csv")), DataFrame(prior_predictive))

posterior_generated_quantities = Chains(generated_quantities(model_generated_quantities, posterior_samples));
mkpath(results_dir("posterior_generated_quantities_csv"))
CSV.write(results_dir("posterior_generated_quantities_csv", savename("posterior_generated_quantities", (@dict max_t), "csv")), DataFrame(posterior_generated_quantities))

Random.seed!(1)
posterior_predictive = predict(model_predict, posterior_samples);
mkpath(results_dir("posterior_predictive_csv"))
CSV.write(results_dir("posterior_predictive_csv", savename("posterior_predictive", (@dict max_t), "csv")), DataFrame(posterior_predictive))