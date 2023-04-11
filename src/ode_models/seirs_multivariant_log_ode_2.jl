# https://q.uiver.app/?q=WzAsMjQsWzMsNCwiU18wIl0sWzcsNCwiSV8xIl0sWzMsOCwiSV8yIl0sWzksNCwiUl8xIl0sWzYsMTAsIlJfMiJdLFs5LDgsIklfezEyfSJdLFszLDYsIkVfMiJdLFs1LDQsIkVfMSJdLFs5LDYsIkVfezEyfSJdLFs2LDE0LCJFX3syMn0iXSxbNCwxNSwiSV97MjJ9Il0sWzcsMiwiSF8xIl0sWzcsMCwiSUNVXzEiXSxbOSwwLCJEXzEiXSxbMTAsMTUsIkhfezIyfSJdLFsxMiwxNSwiSUNVX3syMn0iXSxbMTIsMTcsIkRfezIyfSJdLFsxMCw5LCJIX3sxMn0iXSxbMTEsMTAsIklDVV97MTJ9Il0sWzEyLDExLCJEX3sxMn0iXSxbMiw5LCJIXzIiXSxbMSwxMCwiSUNVXzIiXSxbMCwxMSwiRF8yIl0sWzYsMTIsIlNfezIyfSJdLFsxLDMsIlxcbnVfMSgxLVxcdGV4dHJte2locn1fMSkiXSxbMiw0LCJcXG51XzIoMS1cXHRleHRybXtpaHJ9XzIpIiwyXSxbMCw2LCJcXGJldGFfMihJXzIrSV97MTJ9K0lfezIyfSkiLDJdLFs2LDIsIlxcZ2FtbWFfMiJdLFswLDcsIlxcYmV0YV8xIElfMSJdLFs3LDEsIlxcZ2FtbWFfMSJdLFszLDgsIlxcc2lnbWFfMVxcYmV0YV8yKElfMitJX3sxMn0rSV97MjJ9KSJdLFs4LDUsIlxcZ2FtbWFfMiJdLFs5LDEwLCJcXGdhbW1hXzIiLDJdLFszLDAsIlxca2FwcGFfMSIsMSx7ImN1cnZlIjotNX1dLFs1LDQsIlxcbnVfMigxLVxcdGV4dHJte2locn1fMikiXSxbMSwxMSwiXFxudV8xXFx0ZXh0cm17aWhyfV8xIl0sWzExLDEyLCJcXGV0YV8xXFx0ZXh0cm17aGljdXJ9XzEiXSxbMTIsMTMsIlxcb21lZ2FfMVxcdGV4dHJte2ljdWRyfV8xIl0sWzExLDMsIlxcZXRhXzEoMS1cXHRleHRybXtoaWN1cn1fMSkiXSxbMTIsMywiXFxvbWVnYV8xKDEtXFx0ZXh0cm17aWN1ZHJ9XzEpIl0sWzEwLDE0LCJcXG51XzJcXHRleHRybXtpaHJ9XzMiXSxbMTQsMTUsIlxcZXRhXzJcXHRleHRybXtoaWN1cn1fMyJdLFsxNSwxNiwiXFxvbWVnYV8yXFx0ZXh0cm17aWN1ZHJ9XzMiXSxbNSwxNywiXFxudV8yXFx0ZXh0cm17aWhyfV8yIiwxXSxbMTcsMTgsIlxcZXRhXzJcXHRleHRybXtoaWN1cn1fMiIsMV0sWzE4LDE5LCJcXG9tZWdhXzJcXHRleHRybXtpY3Vkcn1fMiIsMV0sWzE3LDQsIlxcZXRhXzIoMS1cXHRleHRybXtoaWN1cn1fMikiLDFdLFsxOCw0LCJcXG9tZWdhXzIoMS1cXHRleHRybXtpY3Vkcn1fMikiLDFdLFsyLDIwLCJcXG51XzJcXHRleHRybXtpaHJ9XzIiXSxbMjAsMjEsIlxcZXRhXzJcXHRleHRybXtoaWN1cn1fMiJdLFsyMSwyMiwiXFxvbWVnYV8yXFx0ZXh0cm17aWN1ZHJ9XzIiXSxbMjAsNCwiXFxldGFfMigxLVxcdGV4dHJte2hpY3VyfV8yKSJdLFsyMSw0LCJcXG9tZWdhXzIoMS1cXHRleHRybXtpY3Vkcn1fMikiXSxbNCwyMywiXFxrYXBwYV8yIl0sWzIzLDksIlxcYmV0YV8yKElfMitJX3sxMn0rSV97MjJ9KSIsMl0sWzE0LDQsIlxcZXRhXzIoMS1cXHRleHRybXtoaWN1cn1fMykiXSxbMTUsNCwiXFxvbWVnYV8yKDEtXFx0ZXh0cm17aWN1ZHJ9XzMpIl0sWzEwLDQsIlxcbnVfMigxLVxcdGV4dHJte2locn1fMykiLDIseyJjdXJ2ZSI6LTV9XV0=
function seirs_multivariant_log_ode_2!(du, u, p, t)
    (
        S₀, S₂₂,
        E₁, E₁₂, E₂, E₂₂,
        I₁, I₁₂, I₂, I₂₂,
        R₁, R₂,
        H₁, H₁₂, H₂, H₂₂,
        ICU₁, ICU₁₂, ICU₂, ICU₂₂,
        D₁, D₁₂, D₂, D₂₂,
        C₁, C₁₂, C₂, C₂₂
    ) = exp.(u)

    (
        β₁, β₂,
        γ₁, γ₂,
        ν₁, ν₂,
        η₁, η₂,
        ω₁, ω₂,
        κ₁, κ₂,
        ihr₁, ihr₂, ihr₃,
        hicur₁, hicur₂, hicur₃,
        icudr₁, icudr₂, icudr₃,
        σ₁
    ) = p

    N = S₀ + S₂₂ + E₁ + I₁ + R₁ + H₁ + ICU₁ + D₁ + E₂ + I₂ + R₂ + H₂ + ICU₂ + D₂ + E₁₂ + I₁₂ + H₁₂ + ICU₁₂ + D₁₂ + E₂₂ + I₂₂ + H₂₂ + ICU₂₂ + D₂₂

    # Infection
    # S₀ -> E₁
    infection₁ = β₁ * (I₁) * S₀ / N
    # S₀ -> E₂
    infection₂ = β₂ * (I₂ + I₁₂ + I₂₂) * S₀ / N
    # R₁ -> E₁₂
    infection₁₂ = σ₁ * β₂ * (I₂ + I₁₂ + I₂₂) * R₁ / N
    # S₂₂ -> E₂₂
    infection₂₂ = β₂ * (I₂ + I₁₂ + I₂₂) * S₂₂ / N

    # Progression
    # E₁ -> I₁
    progression₁ = γ₁ * E₁
    # E₂ -> I₂
    progression₂ = γ₂ * E₂
    # E₁₂ -> I₁₂
    progression₁₂ = γ₂ * E₁₂
    # E₂₂ -> I₂₂
    progression₂₂ = γ₂ * E₂₂

    # Hospitalization
    # I₁ -> H₁
    hospitalization₁ = ν₁ * ihr₁ * I₁
    # I₂ -> H₂
    hospitalization₂ = ν₂ * ihr₂ * I₂
    # I₁₂ -> H₁₂
    hospitalization₁₂ = ν₂ * ihr₂ * I₁₂
    # I₂₂ -> H₂₂
    hospitalization₂₂ = ν₂ * ihr₃ * I₂₂

    # ICU Admission
    # H₁ -> ICU₁
    icu_admission₁ = η₁ * hicur₁ * H₁
    # H₂ -> ICU₂
    icu_admission₂ = η₂ * hicur₂ * H₂
    # H₁₂ -> ICU₁₂
    icu_admission₁₂ = η₂ * hicur₂ * H₁₂
    # H₂₂ -> ICU₂₂
    icu_admission₂₂ = η₂ * hicur₃ * H₂₂

    # ICU Death
    # ICU₁ -> D₁
    ICU_death₁ = ω₁ * icudr₁ * ICU₁
    # ICU₂ -> D₂
    ICU_death₂ = ω₂ * icudr₂ * ICU₂
    # ICU₁₂ -> D₁₂
    ICU_death₁₂ = ω₂ * icudr₂ * ICU₁₂
    # ICU₂₂ -> D₂₂
    ICU_death₂₂ = ω₂ * icudr₃ * ICU₂₂

    # Non-Hospitalized Recovery
    # I₁ -> R₁
    non_hospitalized_recovery₁ = ν₁ * (1 - ihr₁) * I₁
    # I₂ -> R₂
    non_hospitalized_recovery₂ = ν₂ * (1 - ihr₂) * I₂
    # I₁₂ -> R₂
    non_hospitalized_recovery₁₂ = ν₂ * (1 - ihr₂) * I₁₂
    # I₂₂ -> R₂
    non_hospitalized_recovery₂₂ = ν₂ * (1 - ihr₃) * I₂₂

    # Hospitalized Recovery
    # H₁ -> R₁
    hospitalized_recovery₁ = η₁ * (1 - hicur₁) * H₁
    # H₂ -> R₂
    hospitalized_recovery₂ = η₂ * (1 - hicur₂) * H₂
    # H₁₂ -> R₂
    hospitalized_recovery₁₂ = η₂ * (1 - hicur₂) * H₁₂
    # H₂₂ -> R₂
    hospitalized_recovery₂₂ = η₂ * (1 - hicur₃) * H₂₂

    # ICU Recovery
    # ICU₁ -> R₁
    ICU_recovery₁ = ω₁ * (1 - icudr₁) * ICU₁
    # ICU₂ -> R₂
    ICU_recovery₂ = ω₂ * (1 - icudr₂) * ICU₂
    # ICU₁₂ -> R₂
    ICU_recovery₁₂ = ω₂ * (1 - icudr₂) * ICU₁₂
    # ICU₂₂ -> R₂
    ICU_recovery₂₂ = ω₂ * (1 - icudr₃) * ICU₂₂

    # Waning
    # R₁ -> S₀
    waning₁ = κ₁ * R₁
    # R₂ -> S₂₂
    waning₂ = κ₂ * R₂

    @inbounds begin
        du[1] = (waning₁ - (infection₁ + infection₂)) / S₀
        du[2] = (waning₂ - infection₂₂) / S₂₂
        du[3] = (infection₁ - progression₁) / E₁
        du[4] = (infection₁₂ - progression₁₂) / E₁₂
        du[5] = (infection₂ - progression₂) / E₂
        du[6] = (infection₂₂ - progression₂₂) / E₂₂
        du[7] = (progression₁ - (hospitalization₁ + non_hospitalized_recovery₁)) / I₁
        du[8] = (progression₁₂ - (hospitalization₁₂ + non_hospitalized_recovery₁₂)) / I₁₂
        du[9] = (progression₂ - (hospitalization₂ + non_hospitalized_recovery₂)) / I₂
        du[10] = (progression₂₂ - (hospitalization₂₂ + non_hospitalized_recovery₂₂)) / I₂₂
        du[11] = ((non_hospitalized_recovery₁ + hospitalized_recovery₁ + ICU_recovery₁) - waning₁) / R₁
        du[12] = (((non_hospitalized_recovery₁₂ + hospitalized_recovery₁₂ + ICU_recovery₁₂) +
                   (non_hospitalized_recovery₂ + hospitalized_recovery₂ + ICU_recovery₂) +
                   (non_hospitalized_recovery₂₂ + hospitalized_recovery₂₂ + ICU_recovery₂₂))
                   - waning₂) / R₂
        du[13] = (hospitalization₁ - (icu_admission₁ + hospitalized_recovery₁)) / H₁
        du[14] = (hospitalization₁₂ - (icu_admission₁₂ + hospitalized_recovery₁₂)) / H₁₂
        du[15] = (hospitalization₂ - (icu_admission₂ + hospitalized_recovery₂)) / H₂
        du[16] = (hospitalization₂ - (icu_admission₂₂ + hospitalized_recovery₂₂)) / H₂₂
        du[17] = (icu_admission₁ - (ICU_death₁ + ICU_recovery₁)) / ICU₁
        du[18] = (icu_admission₁₂ - (ICU_death₁₂ + ICU_recovery₁₂)) / ICU₁₂
        du[19] = (icu_admission₂ - (ICU_death₂ + ICU_recovery₂)) / ICU₂
        du[20] = (icu_admission₂₂ - (ICU_death₂₂ + ICU_recovery₂₂)) / ICU₂₂
        du[21] = ICU_death₁ / D₁
        du[22] = ICU_death₁₂ / D₁₂
        du[23] = ICU_death₂ / D₂
        du[24] = ICU_death₂₂ / D₂₂
        du[25] = progression₁ / C₁
        du[26] = progression₁₂ / C₁₂
        du[27] = progression₂ / C₂
        du[28] = progression₂₂ / C₂₂
    end
    nothing
end