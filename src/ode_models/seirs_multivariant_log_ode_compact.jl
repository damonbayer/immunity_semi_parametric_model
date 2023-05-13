# https://q.uiver.app/?q=WzAsMTMsWzAsMCwiU18wIl0sWzEsMCwiRV8xIl0sWzIsMCwiSV8xIl0sWzIsMSwiUl8xIl0sWzEsMiwiRV97Mn0iXSxbMiwyLCJJXzIiXSxbMiwzLCJSXzIiXSxbMCwyLCJTXzIiXSxbMywwLCJIXzEiXSxbNCwwLCJJQ1VfMSJdLFszLDIsIkhfMiJdLFs0LDIsIklDVV8yIl0sWzQsMSwiRCJdLFswLDFdLFsxLDJdLFsyLDNdLFswLDRdLFs0LDVdLFszLDBdLFszLDRdLFs1LDZdLFs2LDddLFsyLDhdLFs4LDldLFs5LDNdLFs4LDNdLFs1LDEwXSxbMTAsNl0sWzEwLDExXSxbMTEsNl0sWzcsNF0sWzExLDEyXSxbOSwxMl1d
function seirs_multivariant_log_ode_compact!(du, u, p, t)
    (
        S₀, S₂,
        E₁, E₂,
        I₁, I₂,
        R₁, R₂,
        H₁, H₂,
        ICU₁, ICU₂,
        D,
        C₁, C₂
    ) = exp.(u)

    (
        β₁, β₂,
        γ₁, γ₂,
        ν₁, ν₂,
        η₁, η₂,
        ω₁, ω₂,
        κ₁, κ₂,
        ihr₁, ihr₂,
        hicur₁, hicur₂,
        icudr₁, icudr₂,
        λ₁
    ) = p

    N = S₀ + S₂ + E₁ + E₂ + I₁ + I₂ + R₁ + R₂ + H₁ + H₂ + ICU₁ + ICU₂ + D

    # Infection
    # S₀ -> E₁
    infection₀₁ = β₁ * I₁ * S₀ / N
    # S₀ -> E₂
    infection₀₂ = β₂ * I₂ * S₀ / N
    # R₁ -> E₂
    infection₁₂ = λ₁ * β₂ * I₂ * R₁ / N
    # S₂ -> E₂
    infection₂₂ = β₂ * I₂ * S₂ / N

    # Progression
    # E₁ -> I₁
    progression₁ = γ₁ * E₁
    # E₂ -> I₂
    progression₂ = γ₂ * E₂

    # Hospitalization
    # I₁ -> H₁
    hospitalization₁ = ν₁ * ihr₁ * I₁
    # I₂ -> H₂
    hospitalization₂ = ν₂ * ihr₂ * I₂

    # ICU Admission
    # H₁ -> ICU₁
    icu_admission₁ = η₁ * hicur₁ * H₁
    # H₂ -> ICU₂
    icu_admission₂ = η₂ * hicur₂ * H₂

    # ICU Death
    # ICU₁ -> D
    ICU_death₁ = ω₁ * icudr₁ * ICU₁
    # ICU₂ -> D
    ICU_death₂ = ω₂ * icudr₂ * ICU₂

    # Non-Hospitalized Recovery
    # I₁ -> R₁
    non_hospitalized_recovery₁ = ν₁ * (1 - ihr₁) * I₁
    # I₂ -> R₂
    non_hospitalized_recovery₂ = ν₂ * (1 - ihr₂) * I₂

    # Hospitalized Recovery
    # H₁ -> R₁
    hospitalized_recovery₁ = η₁ * (1 - hicur₁) * H₁
    # H₂ -> R₂
    hospitalized_recovery₂ = η₂ * (1 - hicur₂) * H₂

    # ICU Recovery
    # ICU₁ -> R₁
    ICU_recovery₁ = ω₁ * (1 - icudr₁) * ICU₁
    # ICU₂ -> R₂
    ICU_recovery₂ = ω₂ * (1 - icudr₂) * ICU₂

    # Waning
    # R₁ -> S₀
    waning₁ = κ₁ * R₁
    # R₂ -> S₂
    waning₂ = κ₂ * R₂

    @inbounds begin
        du[1] =  (waning₁ - (infection₀₁ + infection₀₂)) / S₀
        du[2] =  (waning₂ - infection₂₂) / S₂
        du[3] =  (infection₀₁ - progression₁) / E₁
        du[4] =  (infection₀₂ + infection₁₂ + infection₂₂ - progression₂) / E₂
        du[5] =  (progression₁ - (non_hospitalized_recovery₁ + hospitalization₁)) / I₁
        du[6] =  (progression₂ - (non_hospitalized_recovery₂ + hospitalization₂)) / I₂
        du[7] =  (non_hospitalized_recovery₁ + hospitalized_recovery₁ + ICU_recovery₁ - (waning₁ + infection₁₂)) / R₁
        du[8] =  (non_hospitalized_recovery₂ + hospitalized_recovery₂ + ICU_recovery₂ - waning₂) / R₂
        du[9] =  (hospitalization₁ - (hospitalized_recovery₁ + icu_admission₁)) / H₁
        du[10] = (hospitalization₂ - (hospitalized_recovery₂ + icu_admission₂)) / H₂
        du[11] = (icu_admission₁ - (ICU_recovery₁ + ICU_death₁)) / ICU₁
        du[12] = (icu_admission₂ - (ICU_recovery₂ + ICU_death₂)) / ICU₂
        du[13] = (ICU_death₁ + ICU_death₂) / D
        du[14] = (progression₁) / C₁
        du[15] = (progression₂) / C₂
    end
    nothing
end