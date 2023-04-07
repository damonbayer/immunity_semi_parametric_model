const time_interval_in_days = 7.0

const R₀_init_non_centered_loc = log(1.5)
const dur_latent_non_centered_loc = log(2 / 7)
const dur_infectious_non_centered_loc = log(5 / 7)
const dur_hospitalized_non_centered_loc = log(3 / 7)
const dur_icu_non_centered_loc = log(5 / 7)
const dur_waning_non_centered_loc = log(16)
const IHR_non_centered_loc = logit(0.02)
const HICUR_non_centered_loc = logit(0.15)
const ICUDR_non_centered_loc = logit(0.15)
const E_init_prop_non_centered_loc = logit(700 / 3_000_000)
const I_init_prop_non_centered_loc = logit(1500 / 3_000_000)
const R_init_prop_non_centered_loc = logit(600_000 / 3_000_000)
const ρ_deaths_non_centered_loc = logit(0.9)
const ρ_cases_non_centered_loc = logit(0.2)

const R₀_init_non_centered_scale = 0.2
const dur_latent_non_centered_scale = 0.2
const dur_infectious_non_centered_scale = 0.2
const dur_hospitalized_non_centered_scale = 0.2
const dur_icu_non_centered_scale = 0.2
const dur_waning_non_centered_scale = 0.2
const IHR_non_centered_scale = 0.2
const HICUR_non_centered_scale = 0.2
const ICUDR_non_centered_scale = 0.2
const E_init_prop_non_centered_scale = 0.2
const I_init_prop_non_centered_scale = 0.2
const R_init_prop_non_centered_scale = 0.4
const dur_waning_shape_non_centered_scale = 0.5
const ρ_deaths_non_centered_scale = 0.2
const ρ_cases_non_centered_scale = 0.2

const σ_R₀_non_centered_loc = -2
const σ_R₀_non_centered_scale = 0.1