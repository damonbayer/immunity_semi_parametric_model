if !@isdefined(max_t)
    error("Must define max_t before loading data")
end

if !@isdefined(n_forecast_weeks)
    @warn("n_forecast_weeks is not defined. Assigning n_forecast_weeks = 0")
    n_forecast_weeks = 0
end

if !@isdefined(county_id)
    @warn("county_id is not defined. Assigning county_id = 30")
    county_id = 30
end
# n_forecast_weeks = 10
const time_interval_in_days = 7.0

initialization_values = subset(CSV.read(datadir("initialization_values.csv"), DataFrame), :id => ByRow(x -> x == county_id))
cdph_data = subset(CSV.read(datadir("cdph_data.csv"), DataFrame), :id => ByRow(x -> x == county_id))
seq_data = subset(CSV.read(datadir("seq_data.csv"), DataFrame), :id => ByRow(x -> x == county_id))

const county = initialization_values[1, :county]
const popsize, cases_week_0, H_init, ICU_init, D_init, remaining_population_init = float.(map(x -> initialization_values[1, x], [:population, :cases, :H, :ICU, :D, :remaining_population]))
C_init = 1.0
const data_hospitalizations_full, data_icu_full, data_new_cases_full, data_new_deaths_full = map(x -> cdph_data[:, x], [:hospitalized, :icu, :cases, :deaths])
const data_new_seq_variant_1_full, data_new_seq_variant_2_full = map(x -> seq_data[:, x], [:other_count, :lineage_count])

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
seq_obstimes_full = seq_data[:, :time]
seq_obstimes = seq_obstimes_full[seq_obstimes_full .<= l_data]
seq_obstimes_forecast = seq_obstimes_full[seq_obstimes_full .<= l_data_forecast]
l_seq = length(seq_obstimes)
l_seq_forecast = length(seq_obstimes_forecast)

data_new_seq_variant_1 = data_new_seq_variant_1_full[1:l_seq]
data_new_seq_variant_2 = data_new_seq_variant_2_full[1:l_seq]
data_new_seq_variant_1_forecast = data_new_seq_variant_1_full[1:l_seq_forecast]
data_new_seq_variant_2_forecast = data_new_seq_variant_2_full[1:l_seq_forecast]
missing_data_new_seq_variant_1_forecast = Array{Union{Missing, Int64}}(missing, l_seq_forecast)
missing_data_new_seq_variant_2_forecast = Array{Union{Missing, Int64}}(missing, l_seq_forecast)
