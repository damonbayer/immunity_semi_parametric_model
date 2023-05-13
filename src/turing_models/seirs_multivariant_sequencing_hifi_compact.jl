prob = ODEProblem{true}(seirs_multivariant_log_ode_compact!,
    zeros(19),
    (0.0, obstimes[end]),
    ones(15))

@model function seirs_multivariant_sequencing_hifi_compact(data_new_cases, data_hospitalizations, data_icu, data_new_deaths, data_new_seq_variant_1, data_new_seq_variant_2, obstimes, prob, DEAlgorithm, abstol, reltol, variant_2_import_time, seq_obstimes, init_index_to_report=1)
    l_incidence = length(data_new_deaths)
    l_seq = length(data_new_seq_variant_1)

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

    HICUR₁_uncentered ~ Normal()
    HICUR₂_uncentered ~ Normal()

    ICUDR₁_uncentered ~ Normal()
    ICUDR₂_uncentered ~ Normal()

    λ₁_uncentered ~ Normal()

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

    HICUR₁ = logistic(HICUR₁_uncentered * HICUR₁_scale + HICUR₁_loc)
    HICUR₂ = logistic(HICUR₂_uncentered * HICUR₂_scale + HICUR₂_loc)

    ICUDR₁ = logistic(ICUDR₁_uncentered * ICUDR₁_scale + ICUDR₁_loc)
    ICUDR₂ = logistic(ICUDR₂_uncentered * ICUDR₂_scale + ICUDR₂_loc)

    λ₁ = logistic(λ₁_uncentered * λ₁_scale + λ₁_loc)

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
    #     IHR₁, IHR₂,
    #     HICUR₁, HICUR₂,
    #     ICUDR₁, ICUDR₂,
    #     λ₁
    # ]

    p = [
        β₁, 0,
        γ₁, 0,
        ν₁, 0,
        η₁, 0,
        ω₁, 0,
        κ₁, 0,
        IHR₁, 0,
        HICUR₁, 0,
        ICUDR₁, 0,
        λ₁
    ]

    u0_reg_scale = [
        S₀_init, S₂_init,
        E₁_init, E₂_init,
        I₁_init, I₂_init,
        R₁_init, R₂_init,
        H₁_init, H₂_init,
        ICU₁_init, ICU₂_init,
        D_init,
        C₁_init, C₂_init
    ]

    u0 = log.(u0_reg_scale)

    function affect!(integrator)
        integrator.p[2] = β₂
        integrator.p[4] = γ₂
        integrator.p[6] = ν₂
        integrator.p[8] = η₂
        integrator.p[10] = ω₂
        integrator.p[12] = κ₂
        integrator.p[14] = IHR₂
        integrator.p[16] = HICUR₂
        integrator.p[18] = ICUDR₂
        integrator.u[6] = variant_2_import_value
    end

    param_callback = PresetTimeCallback(variant_2_import_time, affect!, save_positions=(false, false))

    obstimes_with_init = vcat(0, obstimes)
    seq_obstimes_with_init = vcat(2 * seq_obstimes[1] - seq_obstimes[2], seq_obstimes)
    all_solve_times = sort(union(obstimes_with_init, seq_obstimes_with_init))

    sol = solve(prob, DEAlgorithm; callback=param_callback, saveat=all_solve_times, save_start=true, verbose=false, abstol=abstol, reltol=reltol, u0=u0, p=p, tspan=(0.0, obstimes[end]))

    # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
    if !SciMLBase.successful_retcode(sol.retcode)
        Turing.@addlogprob! -Inf
        return
    end

    obstimes_index_in_all_solve_times = [findfirst(==(x), sol.t) for x in obstimes_with_init]

    sol_reg_scale_array = exp.(Array(sol)[:, obstimes_index_in_all_solve_times])
    sol_hospitalizations = vec(sum(sol_reg_scale_array[9:10, 2:end], dims=1))
    sol_icu = vec(sum(sol_reg_scale_array[11:12, 2:end], dims=1))
    sol_new_deaths = vec(sol_reg_scale_array[13, 2:end] - sol_reg_scale_array[13, 1:(end-1)])
    sol_new_cases_variant_1 = vec(sol_reg_scale_array[14, 2:end] - sol_reg_scale_array[14, 1:(end-1)])
    sol_new_cases_variant_2 = vec(sol_reg_scale_array[15, 2:end] - sol_reg_scale_array[15, 1:(end-1)])
    sol_new_cases = sol_new_cases_variant_1 + sol_new_cases_variant_2

    hospitalizations_mean = sol_hospitalizations
    icu_mean = sol_icu
    new_deaths_mean = sol_new_deaths .* ρ_deaths
    new_cases_mean = sol_new_cases * ρ_cases

    for i in init_index_to_report:l_incidence
        data_new_cases[i] ~ NegativeBinomial2(new_cases_mean[i], ϕ_new_cases)
        data_hospitalizations[i] ~ NegativeBinomial2(hospitalizations_mean[i], ϕ_hospitalizations)
        data_icu[i] ~ NegativeBinomial2(icu_mean[i], ϕ_icu)
        data_new_deaths[i] ~ NegativeBinomial2(new_deaths_mean[i], ϕ_new_deaths)
    end

    seq_obstimes_index_in_all_solve_times = [findfirst(==(x), sol.t) for x in seq_obstimes_with_init]
    seq_sol_reg_scale_array = exp.(Array(sol)[:, seq_obstimes_index_in_all_solve_times])

    seq_sol_new_cases_variant_1 = vec(seq_sol_reg_scale_array[14, 2:end] - seq_sol_reg_scale_array[14, 1:(end-1)])
    seq_sol_new_cases_variant_2 = vec(seq_sol_reg_scale_array[15, 2:end] - seq_sol_reg_scale_array[15, 1:(end-1)])
    sol_new_cases = seq_sol_new_cases_variant_1 + seq_sol_new_cases_variant_2

    new_seq_variant_1_mean = seq_sol_new_cases_variant_1 * ρ_cases * ρ_seq
    new_seq_variant_2_mean = seq_sol_new_cases_variant_2 * ρ_cases * ρ_seq

    for i in 1:l_seq
        data_new_seq_variant_1[i] ~ NegativeBinomial2(new_seq_variant_1_mean[i], ϕ_new_seq)
        data_new_seq_variant_2[i] ~ NegativeBinomial2(new_seq_variant_2_mean[i], ϕ_new_seq)
    end

    S₀_compartment = sol_reg_scale_array[1, :]
    S₂_compartment = sol_reg_scale_array[2, :]
    E₁_compartment = sol_reg_scale_array[3, :]
    E₂_compartment = sol_reg_scale_array[4, :]
    I₁_compartment = sol_reg_scale_array[5, :]
    I₂_compartment = sol_reg_scale_array[6, :]
    R₁_compartment = sol_reg_scale_array[7, :]
    R₂_compartment = sol_reg_scale_array[8, :]
    H₁_compartment = sol_reg_scale_array[9, :]
    H₂_compartment = sol_reg_scale_array[10, :]
    ICU₁_compartment = sol_reg_scale_array[11, :]
    ICU₂_compartment = sol_reg_scale_array[12, :]
    D_compartment = sol_reg_scale_array[13, :]
    C₁_compartment = sol_reg_scale_array[14, :]
    C₂_compartment = sol_reg_scale_array[15, :]

    seq_I₁_compartment = seq_sol_reg_scale_array[5, :]
    seq_I₂_compartment = seq_sol_reg_scale_array[6, :]

    dur_latent₁_days = dur_latent₁ * time_interval_in_days
    dur_latent₂_days = dur_latent₂ * time_interval_in_days
    dur_infectious₁_days = dur_infectious₁ * time_interval_in_days
    dur_infectious₂_days = dur_infectious₂ * time_interval_in_days
    dur_hospitalized₁_days = dur_hospitalized₁ * time_interval_in_days
    dur_hospitalized₂_days = dur_hospitalized₂ * time_interval_in_days
    dur_icu₁_days = dur_icu₁ * time_interval_in_days
    dur_icu₂_days = dur_icu₂ * time_interval_in_days
    dur_immune₁_days = dur_immune₁ * time_interval_in_days
    dur_immune₂_days = dur_immune₂ * time_interval_in_days

    return (
        S_init=(S₀_compartment+S₂_compartment)[init_index_to_report],
        E_init=(E₁_compartment+E₂_compartment)[init_index_to_report],
        I_init=(I₁_compartment+I₂_compartment)[init_index_to_report],
        R_init=(R₁_compartment+R₂_compartment)[init_index_to_report],
        H_init=(H₁_compartment+H₂_compartment)[init_index_to_report],
        ICU_init=(ICU₁_compartment+ICU₂_compartment)[init_index_to_report],
        D_init=D_compartment[init_index_to_report],
        C_init=(C₁_compartment+C₂_compartment)[init_index_to_report],
        S_compartment=(S₀_compartment+S₂_compartment)[init_index_to_report:end],
        E_compartment=(E₁_compartment+E₂_compartment)[init_index_to_report:end],
        I_compartment=(I₁_compartment+I₂_compartment)[init_index_to_report:end],
        R_compartment=(R₁_compartment+R₂_compartment)[init_index_to_report:end],
        H_compartment=(H₁_compartment+H₂_compartment)[init_index_to_report:end],
        ICU_compartment=(ICU₁_compartment+ICU₂_compartment)[init_index_to_report:end],
        C_compartment=(C₁_compartment+C₂_compartment)[init_index_to_report:end],
        S₀_compartment=S₀_compartment[init_index_to_report:end],
        S₂_compartment=S₂_compartment[init_index_to_report:end],
        E₁_compartment=E₁_compartment[init_index_to_report:end],
        E₂_compartment=E₂_compartment[init_index_to_report:end],
        I₁_compartment=I₁_compartment[init_index_to_report:end],
        I₂_compartment=I₂_compartment[init_index_to_report:end],
        R₁_compartment=R₁_compartment[init_index_to_report:end],
        R₂_compartment=R₂_compartment[init_index_to_report:end],
        H₁_compartment=H₁_compartment[init_index_to_report:end],
        H₂_compartment=H₂_compartment[init_index_to_report:end],
        ICU₁_compartment=ICU₁_compartment[init_index_to_report:end],
        ICU₂_compartment=ICU₂_compartment[init_index_to_report:end],
        D_compartment=D_compartment[init_index_to_report:end],
        C₁_compartment=C₁_compartment[init_index_to_report:end],
        C₂_compartment=C₂_compartment[init_index_to_report:end],
        prop_variant_2=(I₂_compartment ./ (I₁_compartment+I₂_compartment))[init_index_to_report:end],
        seq_prop_variant_2=seq_I₂_compartment ./ (seq_I₁_compartment + seq_I₂_compartment),
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
        HICUR₁=HICUR₁,
        HICUR₂=HICUR₂,
        ICUDR₁=ICUDR₁,
        ICUDR₂=ICUDR₂,
        λ₁=λ₁,
        ρ_deaths=ρ_deaths,
        ρ_cases=ρ_cases,
        ρ_seq=ρ_seq,
        ϕ_new_cases=ϕ_new_cases,
        ϕ_new_seq=ϕ_new_seq,
        ϕ_hospitalizations=ϕ_hospitalizations,
        ϕ_icu=ϕ_icu,
        ϕ_new_deaths=ϕ_new_deaths,
        hospitalizations_mean=hospitalizations_mean[init_index_to_report:end],
        icu_mean=icu_mean[init_index_to_report:end],
        new_deaths_mean=new_deaths_mean[init_index_to_report:end],
        new_cases_mean=new_cases_mean[init_index_to_report:end],
        new_seq_variant_1_mean=new_seq_variant_1_mean,
        new_seq_variant_2_mean=new_seq_variant_2_mean
    )
end
