using DrWatson
using DifferentialEquations
using Turing
using Random
using DataFrames
using CSV
using LogExpFunctions
using immunity_semi_parametric_model

simulationdir(args...) = projectdir("scripts", "simulation", "sensitivity_analysis_cross_immunity_transmissibility", args...)
datadir(args...) = projectdir("data", "simulation", "sensitivity_analysis_cross_immunity_transmissibility", args...)

mkpath(simulationdir())
mkpath(datadir())

obstimes = collect(1:1:(52*2))
max_index_forecast = length(obstimes)


for takeover_speed in ["slow", "medium", "fast"]
    println(takeover_speed)
    include(simulationdir("shared_constants.txt"))
    include(simulationdir("simulate_data", "prior_constants_simulate_data.jl"))
    include(projectdir("src", "ode_models", "seirs_multivariant_log_ode_compact.jl"))
    include(projectdir("src", "turing_models", "seirs_multivariant_sequencing_hifi_compact.jl"))

    if takeover_speed == "slow"
        variant_2_import_time = variant_2_import_time_slow
        λ₁_loc = logit(0.75)
        R₀_variant_2_loc = R₀_variant_1_loc
    elseif takeover_speed == "medium"
        variant_2_import_time = variant_2_import_time_medium
        λ₁_loc = logit(1.0)
        R₀_variant_2_loc = log(1.75)
    elseif takeover_speed == "fast"
        variant_2_import_time = variant_2_import_time_fast
        λ₁_loc = logit(1.0)
        R₀_variant_2_loc = log(2.5)
    end

    println("R₀_variant_2_loc ", R₀_variant_2_loc)
    println("λ₁_loc ", λ₁_loc)
    println("variant_2_import_time ", variant_2_import_time)


    seq_obstimes = (variant_2_import_time-1):1/7:obstimes[end] # sequences are reported daily, starting one week before the first case
    seq_max_index_forecast = length(seq_obstimes)

    data_zero_cases = zeros(Int64, max_index_forecast)
    data_zero_hospitalizations = zeros(Int64, max_index_forecast)
    data_zero_icu = zeros(Int64, max_index_forecast)
    data_zero_deaths = zeros(Int64, max_index_forecast)
    data_zero_seq_variant_1 = zeros(Int64, seq_max_index_forecast)
    data_zero_seq_variant_2 = zeros(Int64, seq_max_index_forecast)

    data_missing_cases = Array{Union{Missing,Int64}}(missing, max_index_forecast)
    data_missing_hospitalizations = Array{Union{Missing,Int64}}(missing, max_index_forecast)
    data_missing_icu = Array{Union{Missing,Int64}}(missing, max_index_forecast)
    data_missing_deaths = Array{Union{Missing,Int64}}(missing, max_index_forecast)
    data_missing_seq_variant_1 = Array{Union{Missing,Int64}}(missing, seq_max_index_forecast)
    data_missing_seq_variant_2 = Array{Union{Missing,Int64}}(missing, seq_max_index_forecast)

    model_sample = seirs_multivariant_sequencing_hifi_compact(data_zero_cases, data_zero_hospitalizations, data_zero_icu, data_zero_deaths, data_zero_seq_variant_1, data_zero_seq_variant_2, obstimes, prob, RadauIIA5(), 1e-11, 1e-8, variant_2_import_time, seq_obstimes, init_index_to_report)
    model_generated_quantities = seirs_multivariant_sequencing_hifi_compact(data_zero_cases, data_zero_hospitalizations, data_zero_icu, data_zero_deaths, data_zero_seq_variant_1, data_zero_seq_variant_2, obstimes, prob, RadauIIA5(), 1e-13, 1e-10, variant_2_import_time, seq_obstimes, init_index_to_report)
    model_predict = seirs_multivariant_sequencing_hifi_compact(data_missing_cases, data_missing_hospitalizations, data_missing_icu, data_missing_deaths, data_missing_seq_variant_1, data_missing_seq_variant_2, obstimes, prob, RadauIIA5(), 1e-11, 1e-8, variant_2_import_time, seq_obstimes, init_index_to_report)

    n_datasets = 200

    Random.seed!(1)
    prior_samples = sample(model_sample, Prior(), 10)
    prior_summary = summarize(prior_samples)

    true_parameters_chains = Chains(transpose(hcat(repeat([repeat([0.0], length(prior_summary.nt.mean))], n_datasets)...)), prior_summary.nt.parameters)

    Random.seed!(1)
    simulated_data = predict(model_predict, true_parameters_chains)
    CSV.write(datadir(savename("simulated_data", (@dict takeover_speed), "csv")), DataFrame(simulated_data))
    wsave(datadir(savename("simulated_data", (@dict takeover_speed), "jld2")), @dict simulated_data)

    Random.seed!(1)
    true_generated_quantities = Chains(generated_quantities(model_generated_quantities, true_parameters_chains))
    CSV.write(datadir(savename("true_generated_quantities", (@dict takeover_speed), "csv")), DataFrame(true_generated_quantities))
    wsave(datadir(savename("true_generated_quantities", (@dict takeover_speed), "jld2")), @dict true_generated_quantities)
end