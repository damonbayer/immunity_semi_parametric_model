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

const time_interval_in_days = 7.0
const popsize = 3_000_000

simulated_data = load(datadir("simulated_data.jld2"))["simulated_data"][sim_id, :, 1];
true_generated_quantities = load(datadir("true_generated_quantities.jld2"))["true_generated_quantities"][sim_id, :, 1];

prop_variant_2_full = vec(vcat(get(true_generated_quantities, :prop_variant_2)[:prop_variant_2]...))
prop_variant_2_fn = linear_interpolation(0:length(prop_variant_2_full)-1, prop_variant_2_full, extrapolation_bc=Interpolations.Flat())

data_hospitalizations_full, data_icu_full, data_new_cases_variant_1_full, data_new_cases_variant_2_full, data_new_deaths_full, data_new_seq_variant_1_full, data_new_seq_variant_2_full = map(x -> Int.(vec(vcat(get(simulated_data, x)[x]...))), [:data_hospitalizations, :data_icu, :data_new_cases_variant_1, :data_new_cases_variant_2, :data_new_deaths, :data_new_seq_variant_1, :data_new_seq_variant_2])
data_new_cases_full = data_new_cases_variant_1_full + data_new_cases_variant_2_full

const ICU_init, H_init, D_init, C_init =  map(x -> only(vcat(get(true_generated_quantities, x)[x]...)), [:ICU_init, :H_init, :D_init, :C_init])
const remaining_population_init = popsize - (ICU_init + H_init + D_init)

l_data = max_t
l_data_forecast = l_data + n_forecast_weeks

obstimes_forecast = float.(1:l_data_forecast)
obstimes = float.(1:l_data)

param_change_times_forecast = obstimes_forecast[1:end-1]
param_change_times = obstimes[1:end-1]

data_new_cases_variant_1 = data_new_cases_variant_1_full[1:l_data]
data_new_cases_variant_2 = data_new_cases_variant_2_full[1:l_data]
data_new_seq_variant_1 = data_new_seq_variant_1_full[1:l_data]
data_new_seq_variant_2 = data_new_seq_variant_2_full[1:l_data]
data_new_cases = data_new_cases_full[1:l_data]
data_hospitalizations = data_hospitalizations_full[1:l_data]
data_icu = data_icu_full[1:l_data]
data_new_deaths = data_new_deaths_full[1:l_data]

data_new_cases_variant_1_forecast = data_new_cases_variant_1_full[1:l_data_forecast]
data_new_cases_variant_2_forecast = data_new_cases_variant_2_full[1:l_data_forecast]
data_new_seq_variant_1_forecast = data_new_seq_variant_1_full[1:l_data_forecast]
data_new_seq_variant_2_forecast = data_new_seq_variant_2_full[1:l_data_forecast]
data_new_cases_forecast = data_new_cases_full[1:l_data_forecast]
data_hospitalizations_forecast = data_hospitalizations_full[1:l_data_forecast]
data_icu_forecast = data_icu_full[1:l_data_forecast]
data_new_deaths_forecast = data_new_deaths_full[1:l_data_forecast]

missing_data_new_cases_variant_1_forecast = Array{Union{Missing, Int64}}(missing, l_data_forecast)
missing_data_new_cases_variant_2_forecast = Array{Union{Missing, Int64}}(missing, l_data_forecast)
missing_data_new_seq_variant_1_forecast = Array{Union{Missing, Int64}}(missing, l_data_forecast)
missing_data_new_seq_variant_2_forecast = Array{Union{Missing, Int64}}(missing, l_data_forecast)
missing_data_new_cases_forecast = Array{Union{Missing, Int64}}(missing, l_data_forecast)
missing_data_hospitalizations_forecast = Array{Union{Missing, Int64}}(missing, l_data_forecast)
missing_data_icu_forecast = Array{Union{Missing, Int64}}(missing, l_data_forecast)
missing_data_new_deaths_forecast = Array{Union{Missing, Int64}}(missing, l_data_forecast)