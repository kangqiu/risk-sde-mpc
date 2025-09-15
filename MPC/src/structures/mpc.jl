abstract type AbstractMPC
end

abstract type SMPC <:AbstractMPC
end

abstract type StochasticMPC <:AbstractMPC
end

struct DeterministicMPC <: AbstractMPC
    Δt::Int                    #sampling time
    𝒩::UnitRange{<:Int}        #timestep range
    weights::Dict
end

mutable struct RiskMeasureMPC <: StochasticMPC 
    Δt::Int                    #sampling time
    𝒩::UnitRange{<:Int}        #timestep range
    weights::Dict
    risk_measure::AbstractRiskMeasure
    model::ScenarioTree
end

