const min_compartment_size = 1e-5

const R₀₁_mean = log(1.5)
const R₀₂_mean = log(1.5)
const dur_infectious₁_mean = log(5 / 7)
const dur_infectious₂_mean = log(5 / 7)
const dur_latent₁_mean = log(2 / 7)
const dur_latent₂_mean = log(2 / 7)
const dur_hospitalization₁_mean = log(3 / 7)
const dur_hospitalization₂_mean = log(3 / 7)
const dur_icu₁_mean = log(5 / 7)
const dur_icu₂_mean = log(5 / 7)
const ihr₁_mean = logit(0.02)
const ihr₂_mean = logit(0.02)
const ihr₃_mean = logit(0.02)
const hicur₁_mean = logit(0.15)
const hicur₂_mean = logit(0.15)
const hicur₃_mean = logit(0.15)
const icudr₁_mean = logit(0.15)
const icudr₂_mean = logit(0.15)
const icudr₃_mean = logit(0.15)
const σ₁_mean = logit(0.8)
const σ₂_mean = logit(0.8)
const μ_mean = log(0.000001)
const dur_waning_mean = log(20)
const ρ_deaths_mean = logit(0.9)
const ρ_cases_mean = logit(0.2)
const ϕ_new_cases_mean = sqrt(1/100)
const ϕ_new_deaths_mean = sqrt(1/100)

const R₀₁_sd = 0.2
const R₀₂_sd = 0.2
const dur_infectious₁_sd = 0.2
const dur_infectious₂_sd = 0.2
const dur_latent₁_sd = 0.2
const dur_latent₂_sd = 0.2
const dur_hospitalization₁_sd = 0.2
const dur_hospitalization₂_sd = 0.2
const dur_icu₁_sd = 0.2
const dur_icu₂_sd = 0.2
const ihr₁_sd = 0.2
const ihr₂_sd = 0.2
const ihr₃_sd = 0.2
const hicur₁_sd = 0.2
const hicur₂_sd = 0.2
const hicur₃_sd = 0.2
const icudr₁_sd = 0.2
const icudr₂_sd = 0.2
const icudr₃_sd = 0.2
const σ₁_sd = 0.2
const σ₂_sd = 0.2
const μ_sd = 0.27
const dur_waning_sd = 0.2
const ρ_deaths_sd = 0.2
const ρ_cases_sd = 0.2
const ϕ_new_cases_sd = 0.2
const ϕ_new_deaths_sd = 0.2


# Initial conditions
const S₀_init = 1_900_000 - 200
const E₁_init = 100
const I₁_init = 300
const R₁_init = 1_100_000 - 200
const H₁_init = min_compartment_size
const ICU₁_init = min_compartment_size
const D₁_init = min_compartment_size
const C₁_init = min_compartment_size
const E₂_init = min_compartment_size
const I₂_init = min_compartment_size
const R₂_init = min_compartment_size
const H₂_init = min_compartment_size
const ICU₂_init = min_compartment_size
const D₂_init = min_compartment_size
const C₂_init = min_compartment_size
const E₁₂_init = min_compartment_size
const I₁₂_init = min_compartment_size
const R₁₂_init = min_compartment_size
const H₁₂_init = min_compartment_size
const ICU₁₂_init = min_compartment_size
const D₁₂_init = min_compartment_size
const C₁₂_init = min_compartment_size
const E₂₁_init = min_compartment_size
const I₂₁_init = min_compartment_size
const R₂₁_init = min_compartment_size
const H₂₁_init = min_compartment_size
const ICU₂₁_init = min_compartment_size
const D₂₁_init = min_compartment_size
const C₂₁_init = min_compartment_size