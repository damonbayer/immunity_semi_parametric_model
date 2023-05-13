if !@isdefined(max_t)
    error("Must define max_t before loading data")
end

if !@isdefined(n_forecast_weeks)
    @warn("n_forecast_weeks is not defined. Assigning n_forecast_weeks = 0")
    n_forecast_weeks = 0
end

if !@isdefined(sim_id)
    @warn("sim_id is not defined. Assigning sim_id = 1")
    sim_id = 1
end

if data_takeover_speed == "slow"
    variant_2_import_time = variant_2_import_time_slow
elseif data_takeover_speed == "medium"
    variant_2_import_time = variant_2_import_time_medium
elseif data_takeover_speed == "fast"
    variant_2_import_time = variant_2_import_time_fast
end

const time_interval_in_days = 7.0
const popsize = 3_000_000

takeover_speed_dict = Dict(:takeover_speed => data_takeover_speed)
simulated_data = load(datadir(savename("simulated_data", takeover_speed_dict, "jld2")))["simulated_data"][sim_id, :, 1];
true_generated_quantities = load(datadir(savename("true_generated_quantities", takeover_speed_dict, "jld2")))["true_generated_quantities"][sim_id, :, 1];

data_hospitalizations_full, data_icu_full, data_new_cases_full, data_new_deaths_full, data_new_seq_variant_1_full, data_new_seq_variant_2_full = map(x -> Int.(vec(vcat(get(simulated_data, x)[x]...))), [:data_hospitalizations, :data_icu, :data_new_cases, :data_new_deaths, :data_new_seq_variant_1, :data_new_seq_variant_2])

const ICU_init, H_init, D_init, C_init =  map(x -> only(vcat(get(true_generated_quantities, x)[x]...)), [:ICU_init, :H_init, :D_init, :C_init])
const remaining_population_init = popsize - (ICU_init + H_init + D_init)

l_data = max_t
l_data_forecast = l_data + n_forecast_weeks

obstimes_forecast = float.(1:l_data_forecast)
obstimes = float.(1:l_data)

param_change_times_forecast = obstimes_forecast[1:end-1]
param_change_times = obstimes[1:end-1]

data_new_cases = data_new_cases_full[1:l_data]
data_hospitalizations = data_hospitalizations_full[1:l_data]
data_icu = data_icu_full[1:l_data]
data_new_deaths = data_new_deaths_full[1:l_data]

data_new_cases_forecast = data_new_cases_full[1:l_data_forecast]
data_hospitalizations_forecast = data_hospitalizations_full[1:l_data_forecast]
data_icu_forecast = data_icu_full[1:l_data_forecast]
data_new_deaths_forecast = data_new_deaths_full[1:l_data_forecast]

missing_data_new_cases_forecast = Array{Union{Missing, Int64}}(missing, l_data_forecast)
missing_data_hospitalizations_forecast = Array{Union{Missing, Int64}}(missing, l_data_forecast)
missing_data_icu_forecast = Array{Union{Missing, Int64}}(missing, l_data_forecast)
missing_data_new_deaths_forecast = Array{Union{Missing, Int64}}(missing, l_data_forecast)


# sequences are reported daily, starting one week before the first case
seq_obstimes = float.((variant_2_import_time - init_index_to_report):1/7:obstimes[end])
seq_obstimes_forecast = float.((variant_2_import_time - init_index_to_report):1/7:obstimes_forecast[end])
l_seq = length(seq_obstimes)
l_seq_forecast = length(seq_obstimes_forecast)

data_new_seq_variant_1 = data_new_seq_variant_1_full[1:l_seq]
data_new_seq_variant_2 = data_new_seq_variant_2_full[1:l_seq]
data_new_seq_variant_1_forecast = data_new_seq_variant_1_full[1:l_seq_forecast]
data_new_seq_variant_2_forecast = data_new_seq_variant_2_full[1:l_seq_forecast]
missing_data_new_seq_variant_1_forecast = Array{Union{Missing, Int64}}(missing, l_seq_forecast)
missing_data_new_seq_variant_2_forecast = Array{Union{Missing, Int64}}(missing, l_seq_forecast)
