# https://q.uiver.app/?q=WzAsMjYsWzAsNiwiU197M30iXSxbMiw1LCJJXzEiXSxbMiw3LCJJXzIiXSxbNyw1LCJSXzEiXSxbNyw3LCJSXzIiXSxbMTAsNSwiSV97MTJ9Il0sWzEwLDcsIklfezIxfSJdLFsxNSw1LCJSX3sxMn0iXSxbMSw3LCJFXzIiXSxbMSw1LCJFXzEiXSxbOSw1LCJFX3sxMn0iXSxbOSw3LCJFX3syMX0iXSxbMCw4LCJYIl0sWzQsOSwiSF8yIl0sWzYsMTEsIklDVV8yIl0sWzcsMTIsIkRfMiJdLFsxMiw5LCJIX3syMX0iXSxbMTQsMTEsIklDVV97MjF9Il0sWzE1LDEyLCJEX3syMX0iXSxbNCwzLCJIXzEiXSxbNiwxLCJJQ1VfMSJdLFs3LDAsIkRfMSJdLFsxMiwzLCJIX3sxMn0iXSxbMTQsMSwiSUNVX3sxMn0iXSxbMTUsMCwiRF97MTJ9Il0sWzE1LDcsIlJfezIxfSJdLFsxLDMsIlxcbnVfMSgxLVxcdGV4dHJte2locn1fMSkiXSxbMiw0LCJcXG51XzIoMS1cXHRleHRybXtpaHJ9XzIpIiwyXSxbMCw4LCJcXGJldGFfMihJXzIrSV97MTJ9KSIsMl0sWzgsMiwiXFxnYW1tYV8yIl0sWzAsOSwiXFxiZXRhXzEoSV8xK0lfezIxfSkiXSxbOSwxLCJcXGdhbW1hXzEiXSxbMywxMCwiXFxzaWdtYV8yXFxiZXRhXzIoSV8yK0lfezEyfSkiXSxbMTAsNSwiXFxnYW1tYV8yIl0sWzQsMTEsIlxcc2lnbWFfMVxcYmV0YV8xKElfMStJX3syMX0pIiwyXSxbMTEsNiwiXFxnYW1tYV8xIiwyXSxbNSw3LCJcXG51XzIoMS1cXHRleHRybXtpaHJ9XzMpIl0sWzEyLDgsIlxcbXUiLDJdLFsyLDEzLCJcXG51XzJcXHRleHRybXtpaHJ9XzIiXSxbMTMsMTQsIlxcZXRhXzJcXHRleHRybXtoaWN1cn1fMiJdLFsxNCwxNSwiXFxvbWVnYV8yXFx0ZXh0cm17aWN1ZHJ9XzIiXSxbMTMsNCwiXFxldGFfMigxLVxcdGV4dHJte2hpY3VyfV8yKSIsMV0sWzE0LDQsIlxcb21lZ2FfMigxLVxcdGV4dHJte2ljdWRyfV8yKSIsMV0sWzYsMTYsIlxcbnVfMVxcdGV4dHJte2locn1fMyIsMV0sWzE2LDE3LCJcXGV0YV8xXFx0ZXh0cm17aGljdXJ9XzMiLDFdLFsxNywxOCwiXFxvbWVnYV8xXFx0ZXh0cm17aWN1ZHJ9XzMiLDFdLFsxLDE5LCJcXG51XzFcXHRleHRybXtpaHJ9XzEiLDJdLFsxOSwyMCwiXFxldGFfMVxcdGV4dHJte2hpY3VyfV8xIiwyXSxbMjAsMjEsIlxcb21lZ2FfMVxcdGV4dHJte2ljdWRyfV8xIiwyXSxbMjAsMywiXFxvbWVnYV8xKDEtXFx0ZXh0cm17aWN1ZHJ9XzEpIiwxXSxbMTksMywiXFxldGFfMSgxLVxcdGV4dHJte2hpY3VyfV8xKSIsMV0sWzUsMjIsIlxcbnVfMlxcdGV4dHJte2locn1fMyIsMV0sWzIyLDIzLCJcXGV0YV8yXFx0ZXh0cm17aGljdXJ9XzMiLDFdLFsyMyw3LCJcXG9tZWdhXzIoMS1cXHRleHRybXtpY3Vkcn1fMykiLDFdLFsyMiw3LCJcXGV0YV8yKDEtXFx0ZXh0cm17aGljdXJ9XzMpIiwxXSxbMjMsMjQsIlxcb21lZ2FfMlxcdGV4dHJte2ljdWRyfV8zIiwxXSxbMTcsMjUsIlxcb21lZ2FfMSgxLVxcdGV4dHJte2ljdWRyfV8zKSIsMV0sWzE2LDI1LCJcXGV0YV8xKDEtXFx0ZXh0cm17aGljdXJ9XzMpIiwxXSxbNiwyNSwiXFxudV8xKDEtXFx0ZXh0cm17aWhyfV8zKSIsMl0sWzI1LDMsIiIsMCx7ImN1cnZlIjotNH1dLFs3LDMsIiIsMCx7ImN1cnZlIjotNH1dLFszLDBdXQ==
function seirs_multivariant_log_ode!(du, u, p, t)
    (
        S₀,
        E₁, I₁, R₁, H₁, ICU₁, D₁, C₁,
        E₂, I₂, R₂, H₂, ICU₂, D₂, C₂,
        E₁₂, I₁₂, R₁₂, H₁₂, ICU₁₂, D₁₂, C₁₂,
        E₂₁, I₂₁, R₂₁, H₂₁, ICU₂₁, D₂₁, C₂₁
        ) = exp.(u)

    (
        β₁, γ₁, ν₁, ihr₁, η₁, hicur₁, ω₁, icudr₁,
        β₂, γ₂, ν₂, ihr₂, η₂, hicur₂, ω₂, icudr₂,
        ihr₃, hicur₃, icudr₃,
        σ₁, σ₂, μ, waning_rate) = p

    N =  S₀ + E₁ + I₁ + R₁ + H₁ + ICU₁ + D₁ + E₂ + I₂ + R₂ + H₂ + ICU₂ + D₂ + E₁₂ + I₁₂ + R₁₂ + H₁₂ + ICU₁₂ + D₁₂ + E₂₁ + I₂₁ + R₂₁ + H₂₁ + ICU₂₁ + D₂₁

    # S₀ -> E₁
    infection₁ = β₁ * (I₁ + I₂₁) * S₀ / N
    # S₀ -> E₂
    infection₂ = β₂ * (I₂ + I₁₂) * S₀ / N
    # X -> E₂
    importation₂ = μ
    # E₁ -> I₁
    progression₁ = γ₁ * E₁
    # E₂ -> I₂
    progression₂ = γ₂ * E₂
    # I₁ -> H₁
    hospitalization₁ = ν₁ * ihr₁ * I₁
    # I₂ -> H₂
    hospitalization₂ = ν₂ * ihr₂ * I₂
    # I₁ -> R₁
    non_hospitalized_recovery₁ = ν₁ * (1 - ihr₁) * I₁
    # I₂ -> R₂
    non_hospitalized_recovery₂ = ν₂ * (1 - ihr₂) * I₂
    # H₁ -> ICU₁
    icu_admission₁ = η₁ * hicur₁ * H₁
    # H₂ -> ICU₂
    icu_admission₂ = η₂ * hicur₂ * H₂
    # H₁ -> R₁
    hospitalized_recovery₁ = η₁ * (1 - hicur₁) * H₁
    # H₂ -> R₂
    hospitalized_recovery₂ = η₂ * (1 - hicur₂) * H₂
    # ICU₁ -> D₁
    ICU_death₁ = ω₁ * icudr₁ * ICU₁
    # ICU₂ -> D₂
    ICU_death₂ = ω₂ * icudr₂ * ICU₂
    # ICU₁ -> R₁
    ICU_recovery₁ = ω₁ * (1 - icudr₁) * ICU₁
    # ICU₂ -> R₂
    ICU_recovery₂ = ω₂ * (1 - icudr₂) * ICU₂
    # R₁ -> E₁₂
    infection₁₂ = σ₂ * β₂ * (I₂ + I₁₂) * R₁ / N
    # R₂ -> E₂₁
    infection₂₁ = σ₁ * β₁ * (I₁ + I₂₁) * R₂ / N
    # E₁₂ -> I₁₂
    progression₁₂ = γ₂ * E₁₂
    # E₂₁ -> I₂₁
    progression₂₁ = γ₁ * E₂₁
    # I₁₂ -> H₁₂
    hospitalization₁₂ = ν₂ * ihr₃ * I₁₂
    # I₂₁ -> H₂₁
    hospitalization₂₁ = ν₁ * ihr₃ * I₂₁
    # I₁₂ -> R₁₂
    non_hospitalized_recovery₁₂ = ν₂ * (1 - ihr₃) * I₁₂
    # I₂₁ -> R₂₁
    non_hospitalized_recovery₂₁ = ν₁ * (1 - ihr₃) * I₂₁
    # H₁₂ -> ICU₁₂
    icu_admission₁₂ = η₂ * hicur₃ * H₁₂
    # H₂₁ -> ICU₂₁
    icu_admission₂₁ = η₁ * hicur₃ * H₂₁
    # H₁₂ -> R₁₂
    hospitalized_recovery₁₂ = η₂ * (1 - hicur₃) * H₁₂
    # H₂₁ -> R₂₁
    hospitalized_recovery₂₁ = η₁ * (1 - hicur₃) * H₂₁
    # ICU₁₂ -> D₁₂
    ICU_death₁₂ = ω₂ * icudr₃ * ICU₁₂
    # ICU₂₁ -> D₂₁
    ICU_death₂₁ = ω₁ * icudr₃ * ICU₂₁
    # ICU₁₂ -> R₁₂
    ICU_recovery₁₂ = ω₂ * (1 - icudr₃) * ICU₁₂
    # ICU₂₁ -> R₂₁
    ICU_recovery₂₁ = ω₁ * (1 - icudr₃) * ICU₂₁
    # R₂₁ -> R₁
    waning₂₁ = waning_rate * R₂₁
    # R₁₂ -> R₁
    waning₁₂ = waning_rate * R₁₂
    # R₂ -> R₁
    waning₂ = waning_rate * R₂
    # R₁ -> S₀
    waning₁ = waning_rate * R₁ / 1.5

    @inbounds begin
        du[1] = (waning₁ - (infection₁ + infection₂)) / S₀
        du[2] = (infection₁ - progression₁) / E₁
        du[3] = (progression₁ - (hospitalization₁ + non_hospitalized_recovery₁)) / I₁
        du[4] = ((waning₂₁ +  waning₁₂ + waning₂ + non_hospitalized_recovery₁ + hospitalized_recovery₁ + ICU_recovery₁) - (infection₁₂ + waning₁)) / R₁
        du[5] = (hospitalization₁ - (hospitalized_recovery₁ + icu_admission₁)) / H₁
        du[6] = (icu_admission₁ - (ICU_recovery₁ + ICU_death₁)) / ICU₁
        du[7] = ICU_death₁ / D₁
        du[8] = progression₁ / C₁
        du[9] = (infection₂ + importation₂ - progression₂) / E₂
        du[10] = (progression₂ - (hospitalization₂ + non_hospitalized_recovery₂)) / I₂
        du[11] = ((non_hospitalized_recovery₂ + hospitalized_recovery₂ + ICU_recovery₂) - (infection₂₁ + waning₂)) / R₂
        du[12] = (hospitalization₂ - (hospitalized_recovery₂ + icu_admission₂)) / H₂
        du[13] = (icu_admission₂ - (ICU_recovery₂ + ICU_death₂)) / ICU₂
        du[14] = ICU_death₂ / D₂
        du[15] = progression₂ / C₂
        du[16] = (infection₁₂ - progression₁₂) / E₁₂
        du[17] = (progression₁₂ - (hospitalization₁₂ + non_hospitalized_recovery₁₂)) / I₁₂
        du[18] = ((non_hospitalized_recovery₁₂ + hospitalized_recovery₁₂ + ICU_recovery₁₂) - waning₁₂) / R₁₂
        du[19] = (hospitalization₁₂ - (hospitalized_recovery₁₂ + icu_admission₁₂)) / H₁₂
        du[20] = (icu_admission₁₂ - (ICU_recovery₁₂ + ICU_death₁₂)) / ICU₁₂
        du[21] = ICU_death₁₂ / D₁₂
        du[22] = progression₁₂ / C₁₂
        du[23] = (infection₂₁ - progression₂₁) / E₂₁
        du[24] = (progression₂₁ - (hospitalization₂₁ + non_hospitalized_recovery₂₁)) / I₂₁
        du[25] = ((non_hospitalized_recovery₂₁ + hospitalized_recovery₂₁ + ICU_recovery₂₁) - waning₂₁) / R₂₁
        du[26] = (hospitalization₂₁ - (hospitalized_recovery₂₁ + icu_admission₂₁)) / H₂₁
        du[27] = (icu_admission₂₁ - (ICU_recovery₂₁ + ICU_death₂₁)) / ICU₂₁
        du[28] = ICU_death₂₁ / D₂₁
        du[29] = progression₂₁ / C₂₁
    end
    nothing
end