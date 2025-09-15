abstract type DataStruct
end

# struct SDEParameters <: DataStruct
#     θ_v :: Real
#     θ_R :: Real
#     θ_P :: Real
#     θ_Q :: Real

#     ρ :: Real
#     μ :: Real

#     γ₁ :: Real
#     γ₂ :: Real


#     ξ₁ :: Real
#     ξ₂ :: Real
#     ξ₃ :: Real

#     σ_v :: Real
#     σ_R :: Real
#     σ_P :: Real
#     σ_Q :: Real
    
#     wp :: Any #must be interpolation object but no idea how to reference it
#     dwp :: Any

#     θ :: Tuple

# end

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