prob = ODEProblem{true}(seirs_ode_log!,
    zeros(8),
    (0.0, obstimes[end]),
    ones(9))

@model function seirs_unknown_prop_variant_2_contact_rate(prob, logistic_growth_time_offset, data_new_cases, data_new_deaths, data_hospitalizations, data_icu, data_new_seq_variant_1, data_new_seq_variant_2, obstimes, param_change_times, DEAlgorithm, abstol, reltol)
    l_obstimes = length(obstimes)

    # Priors
    dur_latent_non_centered ~ Normal()
    dur_infectious_non_centered ~ Normal()
    dur_hospitalized_non_centered ~ Normal()
    dur_icu_non_centered ~ Normal()
    IHR_non_centered ~ Normal()
    HICUR_non_centered ~ Normal()
    ICUDR_non_centered ~ Normal()

    R₀_saturated_non_centered ~ Normal()
    R₀_shape_non_centered ~ Normal()
    prop_R₀_mixed_non_centered ~ Normal()

    dur_saturated_waning_non_centered ~ Normal()
    dur_waning_shape_non_centered ~ Normal()
    prop_dur_mixed_waning_non_centered ~ Normal()

    ρ_cases_non_centered ~ Normal()
    ρ_deaths_non_centered ~ Normal()
    ρ_seq_non_centered ~ Normal()

    ϕ_hospitalizations_non_centered ~ truncated(Normal(), 0, Inf)
    ϕ_new_cases_non_centered ~ truncated(Normal(), 0, Inf)
    ϕ_new_seq_uncentered ~ truncated(Normal(), 0, Inf)
    ϕ_icu_non_centered ~ truncated(Normal(), 0, Inf)
    ϕ_new_deaths_non_centered ~ truncated(Normal(), 0, Inf)

    prop_variant_2_offset_non_centered ~ Normal()
    variant_ratio_at_logistic_growth_offset_time_non_centered ~ Normal()
    time_to_saturation_non_centered ~ Normal()

    # Would rather do these as a DirichletMultinomial
    E_init_prop_non_centered ~ Normal()
    I_init_prop_non_centered ~ Normal()
    R_init_prop_non_centered ~ Normal()

    # Centering tranformations
    dur_latent = exp(dur_latent_non_centered * dur_latent_non_centered_scale + dur_latent_non_centered_loc)
    dur_infectious = exp(dur_infectious_non_centered * dur_infectious_non_centered_scale + dur_infectious_non_centered_loc)
    dur_hospitalized = exp(dur_hospitalized_non_centered * dur_hospitalized_non_centered_scale + dur_hospitalized_non_centered_loc)
    dur_icu = exp(dur_icu_non_centered * dur_icu_non_centered_scale + dur_icu_non_centered_loc)

    dur_saturated_waning = exp(dur_saturated_waning_non_centered * dur_saturated_waning_non_centered_scale + dur_saturated_waning_non_centered_loc)
    prop_dur_mixed_waning = logistic(prop_dur_mixed_waning_non_centered * prop_dur_mixed_waning_non_centered_scale + prop_dur_mixed_waning_non_centered_loc)
    dur_mixed_waning = dur_saturated_waning * prop_dur_mixed_waning
    dur_waning_shape = exp(dur_waning_shape_non_centered * dur_waning_shape_non_centered_scale + dur_waning_shape_non_centered_loc)

    R₀_saturated = exp(R₀_saturated_non_centered * R₀_saturated_non_centered_scale + R₀_saturated_non_centered_loc)
    prop_R₀_mixed = logistic(prop_R₀_mixed_non_centered * prop_R₀_mixed_non_centered_scale + prop_R₀_mixed_non_centered_loc)
    R₀_mixed = R₀_saturated / prop_R₀_mixed
    R₀_shape = exp(R₀_shape_non_centered * R₀_shape_non_centered_scale + R₀_shape_non_centered_loc)

    ρ_cases = logistic(ρ_cases_non_centered * ρ_cases_non_centered_scale + ρ_cases_non_centered_loc)
    ρ_deaths = logistic(ρ_deaths_non_centered * ρ_deaths_non_centered_scale + ρ_deaths_non_centered_loc)
    ρ_seq = logistic(ρ_seq_non_centered * ρ_seq_non_centered_scale + ρ_seq_non_centered_loc)

    IHR = logistic(IHR_non_centered * IHR_non_centered_scale + IHR_non_centered_loc)
    HICUR = logistic(HICUR_non_centered * HICUR_non_centered_scale + HICUR_non_centered_loc)
    ICUDR = logistic(ICUDR_non_centered * ICUDR_non_centered_scale + ICUDR_non_centered_loc)

    ϕ_hospitalizations = ϕ_hospitalizations_non_centered^-2
    ϕ_new_seq = ϕ_new_seq_uncentered^-2
    ϕ_new_cases = ϕ_new_cases_non_centered^-2
    ϕ_icu = ϕ_icu_non_centered^-2
    ϕ_new_deaths = ϕ_new_deaths_non_centered^-2

    E_init_prop = logistic(E_init_prop_non_centered * E_init_prop_non_centered_scale + E_init_prop_non_centered_loc)
    I_init_prop = logistic(I_init_prop_non_centered * I_init_prop_non_centered_scale + I_init_prop_non_centered_loc)
    R_init_prop = logistic(R_init_prop_non_centered * R_init_prop_non_centered_scale + R_init_prop_non_centered_loc)

    prop_variant_2_offset = exp(prop_variant_2_offset_non_centered * prop_variant_2_offset_non_centered_scale + prop_variant_2_offset_non_centered_loc)

    # variant2:variant1
    variant_ratio_at_logistic_growth_offset_time = logistic(variant_ratio_at_logistic_growth_offset_time_non_centered * variant_ratio_at_logistic_growth_offset_time_non_centered_scale + variant_ratio_at_logistic_growth_offset_time_non_centered_loc)
    time_to_saturation = exp(time_to_saturation_non_centered * time_to_saturation_non_centered_scale + time_to_saturation_non_centered_loc)

    # Logistic growth transformation
    logistic_growth_intercept = log(variant_ratio_at_logistic_growth_offset_time)
    logistic_growth_slope = (logit(0.99) - logit(0.01)) / time_to_saturation

    prop_variant_2 = logistic.(logistic_growth_intercept .+ logistic_growth_slope .* (vcat(0, param_change_times) .- logistic_growth_time_offset))
    prop_variant_2_with_offset = logistic.(logistic_growth_intercept .+ logistic_growth_slope .* (vcat(0, param_change_times) .- logistic_growth_time_offset .- prop_variant_2_offset))

    # Natural scale transformation
    γ = 1 / dur_latent
    ν = 1 / dur_infectious
    η = 1 / dur_hospitalized
    ω = 1 / dur_icu

    R₀_α₀ = log(R₀_saturated)
    R₀_α₁ = (log(R₀_mixed) - log(R₀_saturated)) * 4^R₀_shape

    R₀_t = exp.(R₀_α₀ .+ R₀_α₁ .* (prop_variant_2 .* (1 .- prop_variant_2)) .^ R₀_shape)
    β_t = R₀_t .* ν
    β_init = β_t[1]
    β_t_no_init = β_t[2:end]

    dur_waning_α₀ = log(dur_saturated_waning)
    dur_waning_α₁ = (log(dur_mixed_waning) - log(dur_saturated_waning)) * 4^dur_waning_shape

    dur_waning_t = exp.(dur_waning_α₀ .+ dur_waning_α₁ .* (prop_variant_2_with_offset .* (1 .- prop_variant_2_with_offset)) .^ dur_waning_shape)
    κ_t = 1 ./ dur_waning_t
    κ_init = κ_t[1]
    κ_t_no_init = κ_t[2:end]

    # Initial Condititons
    E_init = E_init_prop * popsize
    I_init = I_init_prop * popsize
    R_init = R_init_prop * popsize
    # H_init is loaded as a constant
    # ICU_init is loaded as a constant
    # D_init is loaded as a constant
    # C_init is loaded as a constant
    # During optimization, E_init, I_init, R_init can sometimes get huge, such that S_init is negative
    # Should reconsider constructing this as a stick-breaking or dirichilet or some other strictly bounded prior
    S_init = max(1, remaining_population_init - (E_init + I_init + R_init))

    u0_reg_scale = [S_init, E_init, I_init, H_init, R_init, C_init, ICU_init, D_init]
    u0 = log.(u0_reg_scale)
    p = [β_init, γ, ν, η, IHR, κ_init, HICUR, ω, ICUDR]

    function param_affect_κ!(integrator)
        ind_t = searchsortedfirst(param_change_times, integrator.t) # Find the index of param_change_times that contains the current timestep
        integrator.p[1] = β_t_no_init[ind_t] # Replace β with a new value from β_t_no_init
        integrator.p[6] = κ_t_no_init[ind_t] # Replace κ with a new value from κ_t_no_init
    end
    param_callback = PresetTimeCallback(param_change_times, param_affect_κ!, save_positions=(false, false))

    sol = solve(prob, DEAlgorithm; callback=param_callback, saveat=obstimes, save_start=true, verbose=false, abstol=abstol, reltol=reltol, u0=u0, p=p, tspan=(0.0, obstimes[end]))

    if sol.retcode != :Success
        Turing.@addlogprob! -Inf
        return
    end

    # Likelihood calculations
    sol_reg_scale_array = exp.(Array(sol))

    sol_hospitalizations = sol_reg_scale_array[4, 2:end]
    sol_new_cases = sol_reg_scale_array[6, 2:end] - sol_reg_scale_array[6, 1:(end-1)]
    sol_icu = sol_reg_scale_array[7, 2:end]
    sol_new_deaths = sol_reg_scale_array[8, 2:end] - sol_reg_scale_array[8, 1:(end-1)]

    hospitalizations_mean = sol_hospitalizations
    new_cases_mean = sol_new_cases .* ρ_cases
    new_seq_mean = new_cases_mean .* ρ_seq
    new_seq_variant_1_mean = new_seq_mean .* (1 .- prop_variant_2)
    new_seq_variant_2_mean = new_seq_mean .* prop_variant_2
    icu_mean = sol_icu
    new_deaths_mean = sol_new_deaths .* ρ_deaths

    for i in 1:l_obstimes
        data_hospitalizations[i] ~ NegativeBinomial2(hospitalizations_mean[i], ϕ_hospitalizations)
        data_new_cases[i] ~ NegativeBinomial2(new_cases_mean[i], ϕ_new_cases)
        data_new_seq_variant_1[i] ~ NegativeBinomial2(new_seq_variant_1_mean[i], ϕ_new_seq)
        data_new_seq_variant_2[i] ~ NegativeBinomial2(new_seq_variant_2_mean[i], ϕ_new_seq)
        data_icu[i] ~ NegativeBinomial2(icu_mean[i], ϕ_icu)
        data_new_deaths[i] ~ NegativeBinomial2(new_deaths_mean[i], ϕ_new_deaths)
    end

    # Generated quantities
    S_compartment = sol_reg_scale_array[1, :]
    E_compartment = sol_reg_scale_array[2, :]
    I_compartment = sol_reg_scale_array[3, :]
    H_compartment = sol_reg_scale_array[4, :]
    R_compartment = sol_reg_scale_array[5, :]
    C_compartment = sol_reg_scale_array[6, :]
    ICU_compartment = sol_reg_scale_array[7, :]
    D_compartment = sol_reg_scale_array[8, :]

    dur_latent_days = dur_latent * time_interval_in_days
    dur_infectious_days = dur_infectious * time_interval_in_days
    dur_hospitalized_days = dur_hospitalized * time_interval_in_days
    dur_waning_days_t = dur_waning_t .* time_interval_in_days
    dur_icu_days = dur_icu * time_interval_in_days

    return (
        IHR=IHR,
        HICUR=HICUR,
        ICUDR=ICUDR,
        β_t=β_t,
        R₀_t=R₀_t,
        R₀_saturated=R₀_saturated,
        prop_R₀_mixed=prop_R₀_mixed,
        R₀_mixed=R₀_mixed,
        R₀_shape=R₀_shape,
        R₀_α₀ = R₀_α₀,
        R₀_α₁ = R₀_α₁,
        dur_latent_days=dur_latent_days,
        dur_infectious_days=dur_infectious_days,
        dur_hospitalized_days=dur_hospitalized_days,
        dur_waning_days_t=dur_waning_days_t,
        dur_icu_days=dur_icu_days,
        S_compartment=S_compartment,
        E_compartment=E_compartment,
        I_compartment=I_compartment,
        H_compartment=H_compartment,
        R_compartment=R_compartment,
        C_compartment=C_compartment,
        ICU_compartment=ICU_compartment,
        D_compartment=D_compartment,
        hospitalizations_mean=hospitalizations_mean,
        new_cases_mean=new_cases_mean,
        new_seq_variant_1_mean=new_seq_variant_1_mean,
        new_seq_variant_2_mean=new_seq_variant_2_mean,
        icu_mean=icu_mean,
        new_deaths_mean=new_deaths_mean,
        dur_saturated_waning=dur_saturated_waning,
        prop_dur_mixed_waning=prop_dur_mixed_waning,
        prop_variant_2_offset=prop_variant_2_offset,
        dur_waning_α₀=dur_waning_α₀,
        dur_waning_α₁=dur_waning_α₁,
        dur_waning_shape=dur_waning_shape,
        variant_ratio_at_logistic_growth_offset_time=variant_ratio_at_logistic_growth_offset_time,
        time_to_saturation=time_to_saturation,
        logistic_growth_intercept=logistic_growth_intercept,
        logistic_growth_slope=logistic_growth_slope,
        ρ_cases=ρ_cases,
        ρ_deaths=ρ_deaths,
        ρ_seq=ρ_seq,
        ϕ_hospitalizations=ϕ_hospitalizations,
        ϕ_new_cases=ϕ_new_cases,
        ϕ_new_seq=ϕ_new_seq,
        ϕ_icu=ϕ_icu,
        ϕ_new_deaths=ϕ_new_deaths
    )
end