using DrWatson
using DifferentialEquations
using Turing
using Random
using DataFrames
using CSV
using LogExpFunctions
using immunity_semi_parametric_model

simulationdir(args...) = projectdir("scripts", "simulation", "compare_inference_unknown_variant_proportion_ihr_hifi", args...)
datadir(args...) = projectdir("data", "simulation", "compare_inference_unknown_variant_proportion_ihr_hifi", args...)

mkpath(simulationdir())
mkpath(datadir())

obstimes = collect(1:1:(52*2))
max_index_forecast = length(obstimes)

include(simulationdir("shared_constants.txt"))
include(simulationdir("simulate_data", "prior_constants_simulate_data.jl"))
include(projectdir("src", "ode_models", "seirs_multivariant_log_ode_2.jl"))
include(projectdir("src", "turing_models", "seirs_multivariant_sequencing_hifi.jl"))

seq_obstimes = (variant_2_import_time - 1):1/7:obstimes[end] # sequences are reported daily, starting one week before the first case
seq_max_index_forecast = length(seq_obstimes)

data_zero_cases_variant_1 = zeros(Int64, max_index_forecast)
data_zero_cases_variant_2 = zeros(Int64, max_index_forecast)
data_zero_seq_variant_1 = zeros(Int64, seq_max_index_forecast)
data_zero_seq_variant_2 = zeros(Int64, seq_max_index_forecast)
data_zero_hospitalizations = zeros(Int64, max_index_forecast)
data_zero_icu = zeros(Int64, max_index_forecast)
data_zero_deaths = zeros(Int64, max_index_forecast)

data_missing_cases_variant_1 = Array{Union{Missing,Int64}}(missing, max_index_forecast)
data_missing_cases_variant_2 = Array{Union{Missing,Int64}}(missing, max_index_forecast)
data_missing_seq_variant_1 = Array{Union{Missing,Int64}}(missing, seq_max_index_forecast)
data_missing_seq_variant_2 = Array{Union{Missing,Int64}}(missing, seq_max_index_forecast)
data_missing_hospitalizations = Array{Union{Missing,Int64}}(missing, max_index_forecast)
data_missing_icu = Array{Union{Missing,Int64}}(missing, max_index_forecast)
data_missing_deaths = Array{Union{Missing,Int64}}(missing, max_index_forecast)

model_sample = seirs_multivariant_sequencing_hifi(data_zero_cases_variant_1, data_zero_cases_variant_2, data_zero_hospitalizations, data_zero_icu, data_zero_deaths, data_zero_seq_variant_1, data_zero_seq_variant_2, obstimes, prob, RadauIIA5(), 1e-11, 1e-8, variant_2_import_time, seq_obstimes, init_index_to_report);
model_generated_quantities = seirs_multivariant_sequencing_hifi(data_zero_cases_variant_1, data_zero_cases_variant_2, data_zero_hospitalizations, data_zero_icu, data_zero_deaths, data_zero_seq_variant_1, data_zero_seq_variant_2, obstimes, prob, RadauIIA5(), 1e-13, 1e-10, variant_2_import_time, seq_obstimes, init_index_to_report);
model_predict = seirs_multivariant_sequencing_hifi(data_missing_cases_variant_1, data_missing_cases_variant_2, data_missing_hospitalizations, data_missing_icu, data_missing_deaths, data_missing_seq_variant_1, data_missing_seq_variant_2, obstimes, prob, RadauIIA5(), 1e-11, 1e-8, variant_2_import_time, seq_obstimes, init_index_to_report);

n_datasets = 200

Random.seed!(1)
prior_samples = sample(model_sample, Prior(), 10)
prior_summary = summarize(prior_samples)

true_parameters_chains = Chains(transpose(hcat(repeat([repeat([0.0], length(prior_summary.nt.mean))], n_datasets)...)), prior_summary.nt.parameters);

Random.seed!(1)
simulated_data = predict(model_predict, true_parameters_chains)
CSV.write(datadir("simulated_data.csv"), DataFrame(simulated_data))
wsave(datadir("simulated_data.jld2"), @dict simulated_data)

Random.seed!(1)
true_generated_quantities = Chains(generated_quantities(model_generated_quantities, true_parameters_chains))
CSV.write(datadir("true_generated_quantities.csv"), DataFrame(true_generated_quantities))
wsave(datadir("true_generated_quantities.jld2"), @dict true_generated_quantities)