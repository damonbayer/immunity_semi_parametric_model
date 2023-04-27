prob = ODEProblem{true}(seirs_multivariant_log_ode_2!,
    zeros(22),
    (0.0, obstimes[end]),
    ones(21))

@model function seirs_multivariant_sequencing(data_new_cases_variant_1, data_new_cases_variant_2, data_hospitalizations, data_icu, data_new_deaths, data_new_seq_variant_1, data_new_seq_variant_2, obstimes, prob, DEAlgorithm, abstol, reltol, init_index_to_report=1)
    l_incidence = length(data_new_deaths)

    R₀_variant_1_uncentered ~ Normal()
    R₀_variant_2_uncentered ~ Normal()

    dur_latent₁_uncentered ~ Normal()
    dur_latent₂_uncentered ~ Normal()

    dur_infectious₁_uncentered ~ Normal()
    dur_infectious₂_uncentered ~ Normal()

    dur_hospitalized₁_uncentered ~ Normal()
    dur_hospitalized₂_uncentered ~ Normal()

    dur_icu₁_uncentered ~ Normal()
    dur_icu₂_uncentered ~ Normal()

    dur_immune₁_uncentered ~ Normal()
    dur_immune₂_uncentered ~ Normal()

    IHR₁_uncentered ~ Normal()
    IHR₂_uncentered ~ Normal()
    IHR₃_uncentered ~ Normal()

    HICUR₁_uncentered ~ Normal()
    HICUR₂_uncentered ~ Normal()
    HICUR₃_uncentered ~ Normal()

    ICUDR₁_uncentered ~ Normal()
    ICUDR₂_uncentered ~ Normal()
    ICUDR₃_uncentered ~ Normal()

    σ₁_uncentered ~ Normal()

    ρ_deaths_uncentered ~ Normal()
    ρ_seq_uncentered ~ Normal()
    ρ_cases_uncentered ~ Normal()

    ϕ_new_cases_uncentered ~ truncated(Normal(), 0, Inf)
    ϕ_new_seq_uncentered ~ truncated(Normal(), 0, Inf)
    ϕ_hospitalizations_uncentered ~ truncated(Normal(), 0, Inf)
    ϕ_icu_uncentered ~ truncated(Normal(), 0, Inf)
    ϕ_new_deaths_uncentered ~ truncated(Normal(), 0, Inf)

    # Centering
    R₀_variant_1 = exp(R₀_variant_1_uncentered * R₀_variant_1_scale + R₀_variant_1_loc)
    R₀_variant_2 = exp(R₀_variant_2_uncentered * R₀_variant_2_scale + R₀_variant_2_loc)

    dur_latent₁ = exp(dur_latent₁_uncentered * dur_latent₁_scale + dur_latent₁_loc)
    dur_latent₂ = exp(dur_latent₂_uncentered * dur_latent₂_scale + dur_latent₂_loc)

    dur_infectious₁ = exp(dur_infectious₁_uncentered * dur_infectious₁_scale + dur_infectious₁_loc)
    dur_infectious₂ = exp(dur_infectious₂_uncentered * dur_infectious₂_scale + dur_infectious₂_loc)

    dur_hospitalized₁ = exp(dur_hospitalized₁_uncentered * dur_hospitalized₁_scale + dur_hospitalized₁_loc)
    dur_hospitalized₂ = exp(dur_hospitalized₂_uncentered * dur_hospitalized₂_scale + dur_hospitalized₂_loc)

    dur_icu₁ = exp(dur_icu₁_uncentered * dur_icu₁_scale + dur_icu₁_loc)
    dur_icu₂ = exp(dur_icu₂_uncentered * dur_icu₂_scale + dur_icu₂_loc)

    dur_immune₁ = exp(dur_immune₁_uncentered * dur_immune₁_scale + dur_immune₁_loc)
    dur_immune₂ = exp(dur_immune₂_uncentered * dur_immune₂_scale + dur_immune₂_loc)

    IHR₁ = logistic(IHR₁_uncentered * IHR₁_scale + IHR₁_loc)
    IHR₂ = logistic(IHR₂_uncentered * IHR₂_scale + IHR₂_loc)
    IHR₃ = logistic(IHR₃_uncentered * IHR₃_scale + IHR₃_loc)

    HICUR₁ = logistic(HICUR₁_uncentered * HICUR₁_scale + HICUR₁_loc)
    HICUR₂ = logistic(HICUR₂_uncentered * HICUR₂_scale + HICUR₂_loc)
    HICUR₃ = logistic(HICUR₃_uncentered * HICUR₃_scale + HICUR₃_loc)

    ICUDR₁ = logistic(ICUDR₁_uncentered * ICUDR₁_scale + ICUDR₁_loc)
    ICUDR₂ = logistic(ICUDR₂_uncentered * ICUDR₂_scale + ICUDR₂_loc)
    ICUDR₃ = logistic(ICUDR₃_uncentered * ICUDR₃_scale + ICUDR₃_loc)

    σ₁ = logistic(σ₁_uncentered * σ₁_scale + σ₁_loc)

    ρ_deaths = logistic(ρ_deaths_uncentered * ρ_deaths_scale + ρ_deaths_loc)
    ρ_cases = logistic(ρ_cases_uncentered * ρ_cases_scale + ρ_cases_loc)
    ρ_seq = logistic(ρ_seq_uncentered * ρ_seq_scale + ρ_seq_loc)

    ϕ_new_cases = (ϕ_new_cases_uncentered * ϕ_new_cases_scale + ϕ_new_cases_loc)^-2
    ϕ_new_seq = (ϕ_new_seq_uncentered * ϕ_new_seq_scale + ϕ_new_seq_loc)^-2
    ϕ_hospitalizations = (ϕ_hospitalizations_uncentered * ϕ_hospitalizations_scale + ϕ_hospitalizations_loc)^-2
    ϕ_icu = (ϕ_icu_uncentered * ϕ_icu_scale + ϕ_icu_loc)^-2
    ϕ_new_deaths = (ϕ_new_deaths_uncentered * ϕ_new_deaths_scale + ϕ_new_deaths_loc)^-2

    # Transformations
    ν₁ = 1 / dur_infectious₁
    ν₂ = 1 / dur_infectious₂
    β₁ = R₀_variant_1 * ν₁
    β₂ = R₀_variant_2 * ν₂
    γ₁ = 1 / dur_latent₁
    γ₂ = 1 / dur_latent₂
    η₁ = 1 / dur_hospitalized₁
    η₂ = 1 / dur_hospitalized₂
    ω₁ = 1 / dur_icu₁
    ω₂ = 1 / dur_icu₂
    κ₁ = 1 / dur_immune₁
    κ₂ = 1 / dur_immune₂

    # p = [
    #     β₁, β₂,
    #     γ₁, γ₂,
    #     ν₁, ν₂,
    #     η₁, η₂,
    #     ω₁, ω₂,
    #     κ₁, κ₂,
    #     IHR₁, IHR₂, IHR₃,
    #     HICUR₁, HICUR₂, HICUR₃,
    #     ICUDR₁, ICUDR₂, ICUDR₃,
    #     σ₁
    # ]

    p = [
        β₁, 0,
        γ₁, 0,
        ν₁, 0,
        η₁, 0,
        ω₁, 0,
        κ₁, 0,
        IHR₁, IHR₂, IHR₃,
        HICUR₁, HICUR₂, HICUR₃,
        ICUDR₁, ICUDR₂, ICUDR₃,
        σ₁
    ]

    u0_reg_scale = [
        S₀_init, S₂₂_init,
        E₁_init, E₁₂_init, E₂_init, E₂₂_init,
        I₁_init, I₁₂_init, I₂_init, I₂₂_init,
        R₁_init, R₂_init,
        H₁_init, H₁₂_init, H₂_init, H₂₂_init,
        ICU₁_init, ICU₁₂_init, ICU₂_init, ICU₂₂_init,
        D₁_init, D₁₂_init, D₂_init, D₂₂_init,
        C₁_init, C₁₂_init, C₂_init, C₂₂_init
    ]

    u0 = log.(u0_reg_scale)

    function affect!(integrator)
        integrator.p[2] = β₂
        integrator.p[4] = γ₂
        integrator.p[6] = ν₂
        integrator.p[8] = η₂
        integrator.p[10] = ω₂
        integrator.p[12] = κ₂
        integrator.u[10] = variant_2_import_value
    end

    param_callback = PresetTimeCallback(variant_2_import_time, affect!, save_positions=(false, false))

    sol = solve(prob, DEAlgorithm; callback=param_callback, saveat=obstimes, save_start=true, verbose=false, abstol=abstol, reltol=reltol, u0=u0, p=p, tspan=(0.0, obstimes[end]))

    # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
    if !SciMLBase.successful_retcode(sol.retcode)
        Turing.@addlogprob! -Inf
        return
    end

    sol_reg_scale_array = exp.(Array(sol))
    sol_hospitalizations = vec(sum(sol_reg_scale_array[13:16, 2:end], dims=1))
    sol_icu = vec(sum(sol_reg_scale_array[17:20, 2:end], dims=1))
    sol_new_deaths = vec(sum(sol_reg_scale_array[21:24, 2:end] - sol_reg_scale_array[21:24, 1:(end-1)], dims=1))
    sol_new_cases_variant_1 = vec(sol_reg_scale_array[25, 2:end] - sol_reg_scale_array[25, 1:(end-1)])
    sol_new_cases_variant_2 = vec(sum(sol_reg_scale_array[26:28, 2:end] - sol_reg_scale_array[26:28, 1:(end-1)], dims=1))

    hospitalizations_mean = sol_hospitalizations
    icu_mean = sol_icu
    new_deaths_mean = sol_new_deaths .* ρ_deaths
    new_cases_variant_1_mean = sol_new_cases_variant_1 .* ρ_cases
    new_cases_variant_2_mean = sol_new_cases_variant_2 .* ρ_cases
    new_seq_variant_1_mean = new_cases_variant_1_mean .* ρ_seq
    new_seq_variant_2_mean = new_cases_variant_2_mean .* ρ_seq

    for i in init_index_to_report:l_incidence
        data_new_cases_variant_1[i] ~ NegativeBinomial2(new_cases_variant_1_mean[i], ϕ_new_cases)
        data_new_cases_variant_2[i] ~ NegativeBinomial2(new_cases_variant_2_mean[i], ϕ_new_cases)
        data_new_seq_variant_1[i] ~ NegativeBinomial2(new_seq_variant_1_mean[i], ϕ_new_seq)
        data_new_seq_variant_2[i] ~ NegativeBinomial2(new_seq_variant_2_mean[i], ϕ_new_seq)
        data_hospitalizations[i] ~ NegativeBinomial2(hospitalizations_mean[i], ϕ_hospitalizations)
        data_icu[i] ~ NegativeBinomial2(icu_mean[i], ϕ_icu)
        data_new_deaths[i] ~ NegativeBinomial2(new_deaths_mean[i], ϕ_new_deaths)
    end

    S₀ = sol_reg_scale_array[1, :]
    S₂₂ = sol_reg_scale_array[2, :]
    E₁ = sol_reg_scale_array[3, :]
    E₁₂ = sol_reg_scale_array[4, :]
    E₂ = sol_reg_scale_array[5, :]
    E₂₂ = sol_reg_scale_array[6, :]
    I₁ = sol_reg_scale_array[7, :]
    I₁₂ = sol_reg_scale_array[8, :]
    I₂ = sol_reg_scale_array[9, :]
    I₂₂ = sol_reg_scale_array[10, :]
    R₁ = sol_reg_scale_array[11, :]
    R₂ = sol_reg_scale_array[12, :]
    H₁ = sol_reg_scale_array[13, :]
    H₁₂ = sol_reg_scale_array[14, :]
    H₂ = sol_reg_scale_array[15, :]
    H₂₂ = sol_reg_scale_array[16, :]
    ICU₁ = sol_reg_scale_array[17, :]
    ICU₁₂ = sol_reg_scale_array[18, :]
    ICU₂ = sol_reg_scale_array[19, :]
    ICU₂₂ = sol_reg_scale_array[20, :]
    D₁ = sol_reg_scale_array[21, :]
    D₁₂ = sol_reg_scale_array[22, :]
    D₂ = sol_reg_scale_array[23, :]
    D₂₂ = sol_reg_scale_array[24, :]
    C₁ = sol_reg_scale_array[25, :]
    C₁₂ = sol_reg_scale_array[26, :]
    C₂ = sol_reg_scale_array[27, :]
    C₂₂ = sol_reg_scale_array[28, :]

    dur_latent₁_days=dur_latent₁ * time_interval_in_days
    dur_latent₂_days=dur_latent₂ * time_interval_in_days
    dur_infectious₁_days=dur_infectious₁ * time_interval_in_days
    dur_infectious₂_days=dur_infectious₂ * time_interval_in_days
    dur_hospitalized₁_days=dur_hospitalized₁ * time_interval_in_days
    dur_hospitalized₂_days=dur_hospitalized₂ * time_interval_in_days
    dur_icu₁_days=dur_icu₁ * time_interval_in_days
    dur_icu₂_days=dur_icu₂ * time_interval_in_days
    dur_immune₁_days=dur_immune₁ * time_interval_in_days
    dur_immune₂_days=dur_immune₂ * time_interval_in_days

    return (
        S_init=(S₀+S₂₂)[init_index_to_report],
        E_init=(E₁+E₁₂+E₂+E₂₂)[init_index_to_report],
        I_init=(I₁+I₁₂+I₂+I₂₂)[init_index_to_report],
        R_init=(R₁+R₂)[init_index_to_report],
        H_init=(H₁+H₁₂+H₂+H₂₂)[init_index_to_report],
        ICU_init=(ICU₁+ICU₁₂+ICU₂+ICU₂₂)[init_index_to_report],
        D_init=(D₁+D₁₂+D₂+D₂₂)[init_index_to_report],
        C_init=(C₁+C₁₂+C₂+C₂₂)[init_index_to_report],
        S=(S₀+S₂₂)[init_index_to_report:end],
        E=(E₁+E₁₂+E₂+E₂₂)[init_index_to_report:end],
        I=(I₁+I₁₂+I₂+I₂₂)[init_index_to_report:end],
        R=(R₁+R₂)[init_index_to_report:end],
        H=(H₁+H₁₂+H₂+H₂₂)[init_index_to_report:end],
        ICU=(ICU₁+ICU₁₂+ICU₂+ICU₂₂)[init_index_to_report:end],
        D=(D₁+D₁₂+D₂+D₂₂)[init_index_to_report:end],
        C=(C₁+C₁₂+C₂+C₂₂)[init_index_to_report:end],
        S₀=S₀[init_index_to_report:end],
        S₂₂=S₂₂[init_index_to_report:end],
        E₁=E₁[init_index_to_report:end],
        E₁₂=E₁₂[init_index_to_report:end],
        E₂=E₂[init_index_to_report:end],
        E₂₂=E₂₂[init_index_to_report:end],
        I₁=I₁[init_index_to_report:end],
        I₁₂=I₁₂[init_index_to_report:end],
        I₂=I₂[init_index_to_report:end],
        I₂₂=I₂₂[init_index_to_report:end],
        R₁=R₁[init_index_to_report:end],
        R₂=R₂[init_index_to_report:end],
        H₁=H₁[init_index_to_report:end],
        H₁₂=H₁₂[init_index_to_report:end],
        H₂=H₂[init_index_to_report:end],
        H₂₂=H₂₂[init_index_to_report:end],
        ICU₁=ICU₁[init_index_to_report:end],
        ICU₁₂=ICU₁₂[init_index_to_report:end],
        ICU₂=ICU₂[init_index_to_report:end],
        ICU₂₂=ICU₂₂[init_index_to_report:end],
        D₁=D₁[init_index_to_report:end],
        D₁₂=D₁₂[init_index_to_report:end],
        D₂=D₂[init_index_to_report:end],
        D₂₂=D₂₂[init_index_to_report:end],
        C₁=C₁[init_index_to_report:end],
        C₁₂=C₁₂[init_index_to_report:end],
        C₂=C₂[init_index_to_report:end],
        C₂₂=C₂₂[init_index_to_report:end],
        prop_variant_2=((I₁₂+I₂+I₂₂)./(I₁+I₁₂+I₂+I₂₂))[init_index_to_report:end],
        R₀_variant_1=R₀_variant_1,
        R₀_variant_2=R₀_variant_2,
        dur_latent₁_days=dur_latent₁_days,
        dur_latent₂_days=dur_latent₂_days,
        dur_infectious₁_days=dur_infectious₁_days,
        dur_infectious₂_days=dur_infectious₂_days,
        dur_hospitalized₁_days=dur_hospitalized₁_days,
        dur_hospitalized₂_days=dur_hospitalized₂_days,
        dur_icu₁_days=dur_icu₁_days,
        dur_icu₂_days=dur_icu₂_days,
        dur_immune₁_days=dur_immune₁_days,
        dur_immune₂_days=dur_immune₂_days,
        IHR₁=IHR₁,
        IHR₂=IHR₂,
        IHR₃=IHR₃,
        HICUR₁=HICUR₁,
        HICUR₂=HICUR₂,
        HICUR₃=HICUR₃,
        ICUDR₁=ICUDR₁,
        ICUDR₂=ICUDR₂,
        ICUDR₃=ICUDR₃,
        σ₁=σ₁,
        ρ_deaths=ρ_deaths,
        ρ_cases=ρ_cases,
        ρ_seq=ρ_seq,
        ϕ_new_cases=ϕ_new_cases,
        ϕ_new_seq=ϕ_new_seq,
        ϕ_hospitalizations=ϕ_hospitalizations,
        ϕ_icu=ϕ_icu,
        ϕ_new_deaths=ϕ_new_deaths,
        IHR_variant_1=IHR₁,
        IHR_variant_2=IHR₂
    )
end
