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
using LinearAlgebra
using FillArrays
using immunity_semi_parametric_model

n_forecast_weeks = 10

immunity_model = "seq-informed"
IHR_model = "seq-informed"
R₀_model = "constant"

datadir(args...) = projectdir("data", "real_data", "ba1", args...)
simulationdir(args...) = projectdir("scripts", "real_data", "ba1_immunity_contact_ihr", args...)
results_dir(args...) = projectdir("results", "real_data", "ba1_immunity_contact_ihr", "SEIRS_ihr_seq_R0_constant_immunity_seq", args...)
mkpath(results_dir())

include(datadir("shared_constants.txt"))
include(simulationdir("load_real_data.jl"))
include(simulationdir("SEIRS_ihr_seq_R0_constant_immunity_seq", "prior_constants.jl"))
include(projectdir("src", "ode_models", "seirs_log_ode.jl"))
include(projectdir("src", "turing_models", "seirs_super_model_hifi.jl"))

model_sample = seirs_super_model_hifi(prob, logistic_growth_time_offset, data_new_cases, data_new_deaths, data_hospitalizations, data_icu, data_new_seq_variant_1, data_new_seq_variant_2, immunity_model, IHR_model, R₀_model, obstimes, param_change_times, seq_obstimes, Tsit5(), 1e-11, 1e-8)
model_optimization = seirs_super_model_hifi(prob, logistic_growth_time_offset, data_new_cases, data_new_deaths, data_hospitalizations, data_icu, data_new_seq_variant_1, data_new_seq_variant_2, immunity_model, IHR_model, R₀_model, obstimes, param_change_times, seq_obstimes, Tsit5(), 1e-13, 1e-10)
model_generated_quantities = seirs_super_model_hifi(prob, logistic_growth_time_offset, data_new_cases_forecast, data_new_deaths_forecast, data_hospitalizations_forecast, data_icu_forecast, data_new_seq_variant_1_forecast, data_new_seq_variant_2_forecast, immunity_model, IHR_model, R₀_model, obstimes_forecast, param_change_times_forecast, seq_obstimes_forecast, Tsit5(), 1e-13, 1e-10)
model_predict = seirs_super_model_hifi(prob, logistic_growth_time_offset, missing_data_new_cases_forecast, missing_data_new_deaths_forecast, missing_data_hospitalizations_forecast, missing_data_icu_forecast, missing_data_new_seq_variant_1_forecast, missing_data_new_seq_variant_2_forecast, immunity_model, IHR_model, R₀_model, obstimes_forecast, param_change_times_forecast, seq_obstimes_forecast, Tsit5(), 1e-13, 1e-10)

prior_samples = load(results_dir("prior_samples_jld2", savename("prior_samples", (@dict max_t), "jld2")))["prior_samples"]
posterior_samples = load(results_dir("posterior_samples_jld2", savename("posterior_samples", (@dict max_t), "jld2")))["posterior_samples"]

prior_samples_augmented = augment_chains_with_forecast_samples(prior_samples, model_sample, model_generated_quantities, "randn")
posterior_samples_augmented = augment_chains_with_forecast_samples(posterior_samples, model_sample, model_generated_quantities, "randn")

prior_generated_quantities = Chains(generated_quantities(model_generated_quantities, prior_samples_augmented));
mkpath(results_dir("prior_generated_quantities_csv"))
CSV.write(results_dir("prior_generated_quantities_csv", savename("prior_generated_quantities", (@dict max_t), "csv")), DataFrame(prior_generated_quantities))

Random.seed!(1)
prior_predictive = predict(model_predict, prior_samples_augmented);
mkpath(results_dir("prior_predictive_csv"))
CSV.write(results_dir("prior_predictive_csv", savename("prior_predictive", (@dict max_t), "csv")), DataFrame(prior_predictive))

posterior_generated_quantities = Chains(generated_quantities(model_generated_quantities, posterior_samples_augmented));
mkpath(results_dir("posterior_generated_quantities_csv"))
CSV.write(results_dir("posterior_generated_quantities_csv", savename("posterior_generated_quantities", (@dict max_t), "csv")), DataFrame(posterior_generated_quantities))

Random.seed!(1)
posterior_predictive = predict(model_predict, posterior_samples_augmented);
mkpath(results_dir("posterior_predictive_csv"))
CSV.write(results_dir("posterior_predictive_csv", savename("posterior_predictive", (@dict max_t), "csv")), DataFrame(posterior_predictive))