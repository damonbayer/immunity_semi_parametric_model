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
n_forecast_weeks = 4
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
IHR_model = model_dict["IHR_model"]

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

model_sample = seirs_super_model_hifi_cdr(prob, logistic_growth_time_offset, data_new_cases, data_new_deaths, data_hospitalizations, data_icu, data_new_seq_variant_1, data_new_seq_variant_2, immunity_model, CDR_model, R₀_model, IHR_model, obstimes, param_change_times, seq_obstimes, Tsit5(), 1e-11, 1e-8)
model_optimization = seirs_super_model_hifi_cdr(prob, logistic_growth_time_offset, data_new_cases, data_new_deaths, data_hospitalizations, data_icu, data_new_seq_variant_1, data_new_seq_variant_2, immunity_model, CDR_model, R₀_model, IHR_model, obstimes, param_change_times, seq_obstimes, Tsit5(), 1e-13, 1e-10)
model_generated_quantities = seirs_super_model_hifi_cdr(prob, logistic_growth_time_offset, data_new_cases_forecast, data_new_deaths_forecast, data_hospitalizations_forecast, data_icu_forecast, data_new_seq_variant_1_forecast, data_new_seq_variant_2_forecast, immunity_model, CDR_model, R₀_model, IHR_model, obstimes_forecast, param_change_times_forecast, seq_obstimes_forecast, Tsit5(), 1e-13, 1e-10)
model_predict = seirs_super_model_hifi_cdr(prob, logistic_growth_time_offset, missing_data_new_cases_forecast, missing_data_new_deaths_forecast, missing_data_hospitalizations_forecast, missing_data_icu_forecast, missing_data_new_seq_variant_1_forecast, missing_data_new_seq_variant_2_forecast, immunity_model, CDR_model, R₀_model, IHR_model, obstimes_forecast, param_change_times_forecast, seq_obstimes_forecast, Tsit5(), 1e-13, 1e-10)

prior_samples = load(results_dir("prior_samples_jld2", savename("prior_samples", (@dict fit_id), "jld2")))["prior_samples"]
posterior_samples = load(results_dir("posterior_samples_jld2", savename("posterior_samples", (@dict fit_id), "jld2")))["posterior_samples"]

durations = DataFrame(wall=MCMCChains.wall_duration(posterior_samples), compute = MCMCChains.compute_duration(posterior_samples))
mkpath(results_dir("duration"))
CSV.write(results_dir("duration", savename("duration", (@dict fit_id), "csv")), durations)

prior_samples_augmented = augment_chains_with_forecast_samples(prior_samples, model_sample, model_generated_quantities, "randn")
posterior_samples_augmented = augment_chains_with_forecast_samples(posterior_samples, model_sample, model_generated_quantities, "randn")

prior_generated_quantities = Chains(generated_quantities(model_generated_quantities, prior_samples_augmented));
mkpath(results_dir("prior_generated_quantities_csv"))
CSV.write(results_dir("prior_generated_quantities_csv", savename("prior_generated_quantities", (@dict fit_id), "csv")), DataFrame(prior_generated_quantities))

Random.seed!(1)
prior_predictive = predict(model_predict, prior_samples_augmented);
mkpath(results_dir("prior_predictive_csv"))
CSV.write(results_dir("prior_predictive_csv", savename("prior_predictive", (@dict fit_id), "csv")), DataFrame(prior_predictive))

posterior_generated_quantities = Chains(generated_quantities(model_generated_quantities, posterior_samples_augmented));
mkpath(results_dir("posterior_generated_quantities_csv"))
CSV.write(results_dir("posterior_generated_quantities_csv", savename("posterior_generated_quantities", (@dict fit_id), "csv")), DataFrame(posterior_generated_quantities))

Random.seed!(1)
posterior_predictive = predict(model_predict, posterior_samples_augmented);
mkpath(results_dir("posterior_predictive_csv"))
CSV.write(results_dir("posterior_predictive_csv", savename("posterior_predictive", (@dict fit_id), "csv")), DataFrame(posterior_predictive))