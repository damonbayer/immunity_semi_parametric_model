# Setup
variant_2_import_time = 35
variant_2_import_value = log(1_000)
min_compartment_size = 1e-5

# Locations
R₀₁_loc = log(1.5)
R₀₂_loc = log(1.5)

dur_latent₁_loc = log(2 / 7)
dur_latent₂_loc = log(2 / 7)

dur_infectious₁_loc =  log(5 / 7)
dur_infectious₂_loc =  log(5 / 7)

dur_hospitalization₁_loc = log(3 / 7)
dur_hospitalization₂_loc = log(3 / 7)

dur_icu₁_loc = log(5 / 7)
dur_icu₂_loc = log(5 / 7)

dur_waning₁_loc = log(20)
dur_waning₂_loc = log(20)

ihr₁_loc = logit(0.02)
ihr₂_loc = logit(0.02)
ihr₃_loc = logit(0.02)

hicur₁_loc = logit(0.15)
hicur₂_loc = logit(0.15)
hicur₃_loc = logit(0.15)

icudr₁_loc = logit(0.15)
icudr₂_loc = logit(0.15)
icudr₃_loc = logit(0.15)

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
R₀₁_scale = 0.2
R₀₂_scale = 0.2

dur_latent₁_scale = 0.2
dur_latent₂_scale = 0.2

dur_infectious₁_scale = 0.2
dur_infectious₂_scale = 0.2

dur_hospitalization₁_scale = 0.2
dur_hospitalization₂_scale = 0.2

dur_icu₁_scale = 0.2
dur_icu₂_scale = 0.2

dur_waning₁_scale = 0.2
dur_waning₂_scale = 0.2

ihr₁_scale = 0.2
ihr₂_scale = 0.2
ihr₃_scale = 0.2

hicur₁_scale = 0.2
hicur₂_scale = 0.2
hicur₃_scale = 0.2

icudr₁_scale = 0.2
icudr₂_scale = 0.2
icudr₃_scale = 0.2

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