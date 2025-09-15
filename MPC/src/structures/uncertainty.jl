struct SDEParameters
    θ_v :: Real
    θ_R :: Real
    θ_P :: Real
    θ_Q :: Real

    ρ :: Real
    μ :: Real

    γ₁ :: Real
    γ₂ :: Real


    ξ₁ :: Real
    ξ₂ :: Real
    ξ₃ :: Real

    σ_v :: Real
    σ_R :: Real
    σ_P :: Real
    σ_Q :: Real
    
    wp :: Any #must be interpolation object but no idea how to reference it
    dwp :: Any

    θ :: Tuple

end

function SDEParameters(wp, dwp)
    θ_v =  0.21
    θ_R = 6.16
    θ_P = 2.88
    θ_Q = 0.29

    ρ = 0.15
    μ = 1.19

    γ₁ = 0.900
    γ₂ = 4.69

    ξ₁ = 0.46
    ξ₂ = 9.47
    ξ₃ = 0.99

    σ_v = exp(-4.11)
    σ_R = exp(1.88)
    σ_P = exp(-2.77)
    σ_Q = exp(0.14)
    
    return SDEParameters(
        θ_v, θ_R, θ_P, θ_Q, ρ, μ, γ₁, γ₂, ξ₁, ξ₂, ξ₃,  σ_v, σ_R, σ_P, σ_Q, wp, dwp, 
        (θ_v, θ_R, θ_P, θ_Q, ρ, μ, γ₁, γ₂, ξ₁, ξ₂, ξ₃,  σ_v, σ_R, σ_P, σ_Q, wp, dwp)
    )
end

abstract type UncertaintyModel
end

struct ScenarioTree <: UncertaintyModel
    ω::Int #number of scenarios
    quantiles::Array{Int}
    # type::String
end

function ScenarioTree(ω::Int)
    #equiprobable quantiles by definition
    seg = ω + 1
    percentiles = collect(1:100)
    L = length(percentiles)
    c = L ÷ seg
    Y = Vector{Vector{Real}}(undef, seg)
    idx = 1
    for i ∈ 1:(seg-1)
        Y[i] = percentiles[idx:idx+c-1]
        idx += c 
    end
    Y[end] = percentiles[idx:end]
    quantiles = [Y[i][end] for i in 1:(length(Y)-1)]

    # if type != "dynamic" || type != "static"
    #     throw("Type of uncertainty model must be static or dynamic")
    # end
    return(ScenarioTree(ω, quantiles))
end

abstract type AbstractRiskMeasure
end

struct Expectation <: AbstractRiskMeasure
end

struct Entropic <: AbstractRiskMeasure
    γ::Real
end

struct VaR <: AbstractRiskMeasure
    α::Int
    objective::Dict
end

struct CVaR <: AbstractRiskMeasure
    α::Int
end

abstract type RiskObjective
end

struct BatteryRisk
end

struct BalanceRisk
end