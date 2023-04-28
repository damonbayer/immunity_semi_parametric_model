const time_interval_in_days = 7.0

# Location
if "seq-informed" ∈ [immunity_model, IHR_model, R₀_model]
    # variant2:variant1 ratio
    const variant_ratio_at_logistic_growth_offset_time_non_centered_loc = logit(0.05)
    const time_to_saturation_non_centered_loc = log(12)
end

if immunity_model == "constant"
    const dur_immune_non_centered_loc = log(20)
elseif immunity_model == "gmrf"
    const σ_dur_immune_non_centered_loc = -2
    const dur_immune_init_non_centered_loc = log(16)
elseif immunity_model == "seq-informed"
    const prop_variant_2_immunity_offset_non_centered_loc = log(1)
    const dur_saturated_immune_non_centered_loc = log(20)
    const prop_dur_mixed_immune_non_centered_loc = logit(0.25)
    const dur_immune_shape_non_centered_loc = log(1)
end

if IHR_model == "constant"
    const IHR_non_centered_loc = logit(0.02)
elseif IHR_model == "gmrf"
    const σ_IHR_non_centered_loc = -2
    const IHR_init_non_centered_loc = logit(0.02)
elseif IHR_model == "seq-informed"
    const IHR_non_centered_variant_1_loc = logit(0.02)
    const IHR_non_centered_variant_2_loc = logit(0.004)
end

if R₀_model == "constant"
    const R₀_non_centered_loc = log(1.5)
elseif R₀_model == "gmrf"
    const σ_R₀_non_centered_loc = -2
    const R₀_init_non_centered_loc = log(1.5)
elseif R₀_model == "seq-informed"
    const R₀_saturated_non_centered_loc = log(1.5)
    const prop_R₀_mixed_non_centered_loc = logit(0.5)
    const R₀_shape_non_centered_loc = log(1)
end

const dur_latent_non_centered_loc = log(2 / 7)
const dur_infectious_non_centered_loc = log(5 / 7)
const dur_hospitalized_non_centered_loc = log(3 / 7)
const dur_icu_non_centered_loc = log(5 / 7)
const ρ_cases_non_centered_loc = logit(0.2)
const ρ_deaths_non_centered_loc = logit(0.9)
const ρ_seq_non_centered_loc = logit(0.1)
const HICUR_non_centered_loc = logit(0.15)
const ICUDR_non_centered_loc = logit(0.15)
const E_init_prop_non_centered_loc = logit(700 / 3_000_000)
const I_init_prop_non_centered_loc = logit(1500 / 3_000_000)
const R_init_prop_non_centered_loc = logit(600_000 / 3_000_000)


# Scale
if "seq-informed" ∈ [immunity_model, IHR_model, R₀_model]
    # variant2:variant1 ratio
    const variant_ratio_at_logistic_growth_offset_time_non_centered_scale = 0.5
    const time_to_saturation_non_centered_scale = 0.4
end

if immunity_model == "constant"
    const dur_immune_non_centered_scale = 0.2
elseif immunity_model == "gmrf"
    const σ_dur_immune_non_centered_scale = 0.1
    const dur_immune_init_non_centered_scale = 0.2
elseif immunity_model == "seq-informed"
    const prop_variant_2_immunity_offset_non_centered_scale = 0.8
    const dur_saturated_immune_non_centered_scale = 0.05
    const prop_dur_mixed_immune_non_centered_scale = 0.3
    const dur_immune_shape_non_centered_scale = 0.5
end

if IHR_model == "constant"
    const IHR_non_centered_scale = 0.2
elseif IHR_model == "gmrf"
    const σ_IHR_non_centered_scale = 0.1
    const IHR_init_non_centered_scale = 0.2
elseif IHR_model == "seq-informed"
    const IHR_non_centered_variant_1_scale = 0.2
    const IHR_non_centered_variant_2_scale = 0.2
end

if R₀_model == "constant"
    const R₀_non_centered_scale = 0.2
elseif R₀_model == "gmrf"
    const σ_R₀_non_centered_scale = 0.1
    const R₀_init_non_centered_scale = 0.2
elseif R₀_model == "seq-informed"
    const R₀_saturated_non_centered_scale = 0.2
    const prop_R₀_mixed_non_centered_scale = 0.3
    const R₀_shape_non_centered_scale = 0.3
end

const dur_latent_non_centered_scale = 0.2
const dur_infectious_non_centered_scale = 0.2
const dur_hospitalized_non_centered_scale = 0.2
const dur_icu_non_centered_scale = 0.2
const ρ_cases_non_centered_scale = 0.2
const ρ_deaths_non_centered_scale = 0.2
const ρ_seq_non_centered_scale = 0.2
const HICUR_non_centered_scale = 0.2
const ICUDR_non_centered_scale = 0.2
const E_init_prop_non_centered_scale = 0.2
const I_init_prop_non_centered_scale = 0.2
const R_init_prop_non_centered_scale = 0.4
