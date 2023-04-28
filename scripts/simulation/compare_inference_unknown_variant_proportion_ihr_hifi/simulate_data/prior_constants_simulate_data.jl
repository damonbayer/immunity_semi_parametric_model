# Setup
variant_2_import_value = log(1_000)
min_compartment_size = 1e-5

time_interval_in_days = 7.0

# Locations
R₀_variant_1_loc = log(1.5)
R₀_variant_2_loc = log(1.5)

dur_latent₁_loc = log(2 / 7)
dur_latent₂_loc = log(2 / 7)

dur_infectious₁_loc =  log(5 / 7)
dur_infectious₂_loc =  log(5 / 7)

dur_hospitalized₁_loc = log(3 / 7)
dur_hospitalized₂_loc = log(3 / 7)

dur_icu₁_loc = log(5 / 7)
dur_icu₂_loc = log(5 / 7)

dur_immune₁_loc = log(20)
dur_immune₂_loc = log(20)

IHR₁_loc = logit(0.02)
IHR₂_loc = logit(0.004)
IHR₃_loc = logit(0.004)

HICUR₁_loc = logit(0.15)
HICUR₂_loc = logit(0.15)
HICUR₃_loc = logit(0.15)

ICUDR₁_loc = logit(0.15)
ICUDR₂_loc = logit(0.15)
ICUDR₃_loc = logit(0.15)

σ₁_loc = logit(1.0)

ρ_deaths_loc = logit(0.9)
ρ_cases_loc = logit(0.2)
ρ_seq_loc = logit(0.1)

ϕ_new_cases_loc = sqrt(1/100)
ϕ_new_seq_loc = sqrt(1/100)
ϕ_hospitalizations_loc = sqrt(1/100)
ϕ_icu_loc = sqrt(1/100)
ϕ_new_deaths_loc = sqrt(1/100)

# Scales
R₀_variant_1_scale = 0.2
R₀_variant_2_scale = 0.2

dur_latent₁_scale = 0.2
dur_latent₂_scale = 0.2

dur_infectious₁_scale = 0.2
dur_infectious₂_scale = 0.2

dur_hospitalized₁_scale = 0.2
dur_hospitalized₂_scale = 0.2

dur_icu₁_scale = 0.2
dur_icu₂_scale = 0.2

dur_immune₁_scale = 0.2
dur_immune₂_scale = 0.2

IHR₁_scale = 0.2
IHR₂_scale = 0.2
IHR₃_scale = 0.2

HICUR₁_scale = 0.2
HICUR₂_scale = 0.2
HICUR₃_scale = 0.2

ICUDR₁_scale = 0.2
ICUDR₂_scale = 0.2
ICUDR₃_scale = 0.2

σ₁_scale = 0.2

ρ_deaths_scale = 0.2
ρ_cases_scale = 0.2
ρ_seq_scale = 0.2

ϕ_hospitalizations_scale = 0.2
ϕ_new_cases_scale = 0.2
ϕ_new_seq_scale = 0.2
ϕ_icu_scale = 0.2
ϕ_new_deaths_scale = 0.2

# Initial conditions
E₁_init = 100
I₁_init = 3 * E₁_init
R₁_init = 1_100_000
S₀_init = 3_000_000 - (E₁_init + I₁_init + R₁_init)
E₁₂_init = min_compartment_size
E₂_init = min_compartment_size
E₂₂_init = min_compartment_size
I₁₂_init = min_compartment_size
I₂_init = min_compartment_size
I₂₂_init = min_compartment_size
R₂_init = min_compartment_size
S₂₂_init = min_compartment_size
H₁_init = min_compartment_size
H₁₂_init = min_compartment_size
H₂_init = min_compartment_size
H₂₂_init = min_compartment_size
ICU₁_init = min_compartment_size
ICU₁₂_init = min_compartment_size
ICU₂_init = min_compartment_size
ICU₂₂_init = min_compartment_size
D₁_init = min_compartment_size
D₁₂_init = min_compartment_size
D₂_init = min_compartment_size
D₂₂_init = min_compartment_size
C₁_init = min_compartment_size
C₁₂_init = min_compartment_size
C₂_init = min_compartment_size
C₂₂_init = min_compartment_size