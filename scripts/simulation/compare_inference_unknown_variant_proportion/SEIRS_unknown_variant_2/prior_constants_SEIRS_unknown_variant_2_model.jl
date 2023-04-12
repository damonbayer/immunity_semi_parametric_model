const time_interval_in_days = 7.0

# Location
const dur_latent_non_centered_loc = log(2 / 7)
const dur_infectious_non_centered_loc = log(5 / 7)
const dur_hospitalized_non_centered_loc = log(3 / 7)
const dur_icu_non_centered_loc = log(5 / 7)

const IHR_non_centered_loc = logit(0.02)
const HICUR_non_centered_loc = logit(0.15)
const ICUDR_non_centered_loc = logit(0.15)

const R₀_saturated_non_centered_loc = log(1.5)
const R₀_shape_non_centered_loc = log(1)
const prop_R₀_mixed_non_centered_loc = logit(0.5)

const dur_saturated_waning_non_centered_loc = log(20)
const dur_waning_shape_non_centered_loc = log(1)
const prop_dur_mixed_waning_non_centered_loc = logit(0.5)

const ρ_cases_non_centered_loc = logit(0.2)
const ρ_deaths_non_centered_loc = logit(0.9)
const ρ_seq_non_centered_loc = logit(0.1)

const prop_variant_2_offset_non_centered_loc = log(1)
const variant_ratio_at_logistic_growth_offset_time_non_centered_loc = logit(0.05)
const time_to_saturation_non_centered_loc = log(12)

const E_init_prop_non_centered_loc = logit(700 / 3_000_000)
const I_init_prop_non_centered_loc = logit(1500 / 3_000_000)
const R_init_prop_non_centered_loc = logit(600_000 / 3_000_000)

# Scale
const dur_latent_non_centered_scale = 0.2
const dur_infectious_non_centered_scale = 0.2
const dur_hospitalized_non_centered_scale = 0.2
const dur_icu_non_centered_scale = 0.2

const IHR_non_centered_scale = 0.2
const HICUR_non_centered_scale = 0.2
const ICUDR_non_centered_scale = 0.2

const R₀_saturated_non_centered_scale = 0.2
const R₀_shape_non_centered_scale = 0.3
const prop_R₀_mixed_non_centered_scale = 0.3

const dur_saturated_waning_non_centered_scale = 0.05
const dur_waning_shape_non_centered_scale = 0.3
const prop_dur_mixed_waning_non_centered_scale = 0.3

const ρ_cases_non_centered_scale = 0.2
const ρ_deaths_non_centered_scale = 0.2
const ρ_seq_non_centered_scale = 0.2

const prop_variant_2_offset_non_centered_scale = 0.8
const variant_ratio_at_logistic_growth_offset_time_non_centered_scale = 0.5
const time_to_saturation_non_centered_scale = 0.4

const E_init_prop_non_centered_scale = 0.2
const I_init_prop_non_centered_scale = 0.2
const R_init_prop_non_centered_scale = 0.4