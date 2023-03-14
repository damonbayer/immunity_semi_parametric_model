prob = ODEProblem{true}(seirs_multivariant_log_ode!,
    zeros(22),
    (0.0, obstimes[end]),
    ones(21))

@model function seirs_multivariant(data_new_cases, data_hospitalizations, data_icu, data_new_deaths, obstimes, prob, DEAlgorithm, abstol, reltol, init_index_to_report=1)
    l_incidence = length(data_new_deaths)

    R₀₁_uncentered ~ Normal()
    R₀₂_uncentered ~ Normal()
    dur_infectious₁_uncentered ~ Normal()
    dur_infectious₂_uncentered ~ Normal()
    dur_latent₁_uncentered ~ Normal()
    dur_latent₂_uncentered ~ Normal()
    dur_hospitalization₁_uncentered ~ Normal()
    dur_hospitalization₂_uncentered ~ Normal()
    dur_icu₁_uncentered ~ Normal()
    dur_icu₂_uncentered ~ Normal()
    ihr₁_uncentered ~ Normal()
    ihr₂_uncentered ~ Normal()
    ihr₃_uncentered ~ Normal()
    hicur₁_uncentered ~ Normal()
    hicur₂_uncentered ~ Normal()
    hicur₃_uncentered ~ Normal()
    icudr₁_uncentered ~ Normal()
    icudr₂_uncentered ~ Normal()
    icudr₃_uncentered ~ Normal()
    σ₁_uncentered ~ Normal()
    σ₂_uncentered ~ Normal()
    μ_uncentered ~ Normal()
    dur_waning_uncentered ~ Normal()
    ρ_deaths_uncentered ~ Normal()
    ρ_cases_uncentered ~ Normal()
    ϕ_new_cases_uncentered ~ truncated(Normal(), 0, Inf)
    ϕ_new_deaths_uncentered ~ truncated(Normal(), 0, Inf)

    # Centering
    R₀₁ = exp(R₀₁_uncentered * R₀₁_sd + R₀₁_mean)
    R₀₂ = exp(R₀₂_uncentered * R₀₂_sd + R₀₂_mean)
    dur_infectious₁ = exp(dur_infectious₁_uncentered * dur_infectious₁_sd + dur_infectious₁_mean)
    dur_infectious₂ = exp(dur_infectious₂_uncentered * dur_infectious₂_sd + dur_infectious₂_mean)
    dur_latent₁ = exp(dur_latent₁_uncentered * dur_latent₁_sd + dur_latent₁_mean)
    dur_latent₂ = exp(dur_latent₂_uncentered * dur_latent₂_sd + dur_latent₂_mean)
    dur_hospitalization₁ = exp(dur_hospitalization₁_uncentered * dur_hospitalization₁_sd + dur_hospitalization₁_mean)
    dur_hospitalization₂ = exp(dur_hospitalization₂_uncentered * dur_hospitalization₂_sd + dur_hospitalization₂_mean)
    dur_icu₁ = exp(dur_icu₁_uncentered * dur_icu₁_sd + dur_icu₁_mean)
    dur_icu₂ = exp(dur_icu₂_uncentered * dur_icu₂_sd + dur_icu₂_mean)
    ihr₁ = logistic(ihr₁_uncentered * ihr₁_sd + ihr₁_mean)
    ihr₂ = logistic(ihr₂_uncentered * ihr₂_sd + ihr₂_mean)
    ihr₃ = logistic(ihr₃_uncentered * ihr₃_sd + ihr₃_mean)
    hicur₁ = logistic(hicur₁_uncentered * hicur₁_sd + hicur₁_mean)
    hicur₂ = logistic(hicur₂_uncentered * hicur₂_sd + hicur₂_mean)
    hicur₃ = logistic(hicur₃_uncentered * hicur₃_sd + hicur₃_mean)
    icudr₁ = logistic(icudr₁_uncentered * icudr₁_sd + icudr₁_mean)
    icudr₂ = logistic(icudr₂_uncentered * icudr₂_sd + icudr₂_mean)
    icudr₃ = logistic(icudr₃_uncentered * icudr₃_sd + icudr₃_mean)
    σ₁ = logistic(σ₁_uncentered * σ₁_sd + σ₁_mean)
    σ₂ = logistic(σ₂_uncentered * σ₂_sd + σ₂_mean)
    μ = exp(μ_uncentered * μ_sd + μ_mean)
    dur_waning = exp(dur_waning_uncentered * dur_waning_sd + dur_waning_mean)
    ρ_deaths = logistic(ρ_deaths_uncentered * ρ_deaths_sd + ρ_deaths_mean)
    ρ_cases = logistic(ρ_cases_uncentered * ρ_cases_sd + ρ_cases_mean)
    ϕ_new_cases = (ϕ_new_cases_uncentered * ϕ_new_cases_sd + ϕ_new_cases_mean)^-2
    ϕ_new_deaths = (ϕ_new_deaths_uncentered * ϕ_new_deaths_sd + ϕ_new_deaths_mean)^-2

    # Transformations
    ν₁ = 1 / dur_infectious₁
    ν₂ = 1 / dur_infectious₂
    β₁ = R₀₁ * ν₁
    β₂ = R₀₂ * ν₂
    γ₁ = 1 / dur_latent₁
    γ₂ = 1 / dur_latent₂
    η₁ = 1 / dur_hospitalization₁
    η₂ = 1 / dur_hospitalization₂
    ω₁ = 1 / dur_icu₁
    ω₂ = 1 / dur_icu₂
    waning_rate = 1 / dur_waning

    p = [β₁, γ₁, ν₁, ihr₁, η₁, hicur₁, ω₁, icudr₁,
        β₂, γ₂, ν₂, ihr₂, η₂, hicur₂, ω₂, icudr₂,
        ihr₃, hicur₃, icudr₃,
        σ₁, σ₂, μ, waning_rate]

    u0_reg_scale = [
        S₀_init,
        E₁_init, I₁_init, R₁_init, H₁_init, ICU₁_init, D₁_init, C₁_init,
        E₂_init, I₂_init, R₂_init, H₂_init, ICU₂_init, D₂_init, C₂_init,
        E₁₂_init, I₁₂_init, R₁₂_init, H₁₂_init, ICU₁₂_init, D₁₂_init, C₁₂_init,
        E₂₁_init, I₂₁_init, R₂₁_init, H₂₁_init, ICU₂₁_init, D₂₁_init, C₂₁_init]

    u0 = log.(u0_reg_scale)

    sol = solve(prob, DEAlgorithm; saveat=obstimes, save_start=true, verbose=false, abstol=abstol, reltol=reltol, u0=u0, p=p, tspan=(0.0, obstimes[end]))

    # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
    if sol.retcode != :Success
        Turing.@addlogprob! -Inf
        return
    end

    sol_reg_scale_array = exp.(Array(sol))

    sol_hospitalizations = sol_reg_scale_array[5, 2:end] +
                           sol_reg_scale_array[12, 2:end] +
                           sol_reg_scale_array[19, 2:end] +
                           sol_reg_scale_array[26, 2:end]

    sol_icu = sol_reg_scale_array[6, 2:end] +
              sol_reg_scale_array[13, 2:end] +
              sol_reg_scale_array[20, 2:end] +
              sol_reg_scale_array[27, 2:end]

    sol_new_deaths =
        sol_reg_scale_array[7, 2:end] - sol_reg_scale_array[7, 1:(end-1)] +
        sol_reg_scale_array[14, 2:end] - sol_reg_scale_array[14, 1:(end-1)] +
        sol_reg_scale_array[21, 2:end] - sol_reg_scale_array[21, 1:(end-1)] +
        sol_reg_scale_array[28, 2:end] - sol_reg_scale_array[28, 1:(end-1)]

    sol_new_cases =
        sol_reg_scale_array[8, 2:end] - sol_reg_scale_array[8, 1:(end-1)] +
        sol_reg_scale_array[15, 2:end] - sol_reg_scale_array[15, 1:(end-1)] +
        sol_reg_scale_array[22, 2:end] - sol_reg_scale_array[22, 1:(end-1)] +
        sol_reg_scale_array[29, 2:end] - sol_reg_scale_array[29, 1:(end-1)]

    hospitalizations_mean = sol_hospitalizations
    icu_mean = sol_icu
    new_deaths_mean = sol_new_deaths .* ρ_deaths
    new_cases_mean = sol_new_cases .* ρ_cases


    sol_reg_scale_array
    for i in init_index_to_report:l_incidence
        data_new_cases[i] ~ NegativeBinomial2(new_cases_mean[i], ϕ_new_cases)
        data_hospitalizations[i] ~ Poisson(hospitalizations_mean[i])
        data_icu[i] ~ Poisson(icu_mean[i])
        data_new_deaths[i] ~ NegativeBinomial2(new_deaths_mean[i], ϕ_new_deaths)
    end

    S₀ = sol_reg_scale_array[1, :]
    E₁ = sol_reg_scale_array[2, :]
    I₁ = sol_reg_scale_array[3, :]
    R₁ = sol_reg_scale_array[4, :]
    H₁ = sol_reg_scale_array[5, :]
    ICU₁ = sol_reg_scale_array[6, :]
    D₁ = sol_reg_scale_array[7, :]
    C₁ = sol_reg_scale_array[8, :]
    E₂ = sol_reg_scale_array[9, :]
    I₂ = sol_reg_scale_array[10, :]
    R₂ = sol_reg_scale_array[11, :]
    H₂ = sol_reg_scale_array[12, :]
    ICU₂ = sol_reg_scale_array[13, :]
    D₂ = sol_reg_scale_array[14, :]
    C₂ = sol_reg_scale_array[15, :]
    E₁₂ = sol_reg_scale_array[16, :]
    I₁₂ = sol_reg_scale_array[17, :]
    R₁₂ = sol_reg_scale_array[18, :]
    H₁₂ = sol_reg_scale_array[19, :]
    ICU₁₂ = sol_reg_scale_array[20, :]
    D₁₂ = sol_reg_scale_array[21, :]
    C₁₂ = sol_reg_scale_array[22, :]
    E₂₁ = sol_reg_scale_array[23, :]
    I₂₁ = sol_reg_scale_array[24, :]
    R₂₁ = sol_reg_scale_array[25, :]
    H₂₁ = sol_reg_scale_array[26, :]
    ICU₂₁ = sol_reg_scale_array[27, :]
    D₂₁ = sol_reg_scale_array[28, :]
    C₂₁ = sol_reg_scale_array[29, :]

    return (
        S_init = (S₀ + R₁)[init_index_to_report],
        E_init = (E₁ + E₂ + E₁₂ + E₂₁)[init_index_to_report],
        I_init = (I₁ + I₂ + I₁₂ + I₂₁)[init_index_to_report],
        R_init = (R₂ + R₁₂ + R₂₁)[init_index_to_report],
        ICU_init = (ICU₁ + ICU₂ + ICU₁₂ + ICU₂₁)[init_index_to_report],
        H_init = (H₁ + H₂ + H₁₂ + H₂₁)[init_index_to_report],
        D_init = (D₁ + D₂ + D₁₂ + D₂₁)[init_index_to_report],
        C_init = (C₁ + C₂ + C₁₂ + C₂₁)[init_index_to_report],
        I₁=I₁[init_index_to_report:end],
        I₂=I₂[init_index_to_report:end],
        I₁₂=I₁₂[init_index_to_report:end],
        I₂₁=I₂₁[init_index_to_report:end],
        prop_variant_2=((I₂ + I₁₂) ./ (I₁ + I₂ + I₁₂ + I₂₁))[init_index_to_report:end])
end
