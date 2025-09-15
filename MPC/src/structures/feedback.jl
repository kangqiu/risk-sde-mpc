abstract type Action end

mutable struct GTAction <: Action 
    const id::Any
    load::Real
end

mutable struct BatteryAction <: Action 
    const id::Any
    charge::Real
    discharge::Real
end

mutable struct WindAction <: Action
    const id::Any
    curtailment::Real
end

mutable struct DemandAction <: Action
    const id::Any
    P::Real
end

mutable struct SimpleMarketAction <: Action
    const id::Any
    buy::Real
    sell::Real
end
mutable struct SpotAction <: Action
    const id::Any
    buy::Vector{<:Real}
    sell::Vector{<:Real}
end


abstract type Feedback
end

mutable struct GTFeedback <: Feedback
    const id::Any
    η::Real
    power::Real
    co2::Real
end

mutable struct BatteryFeedback <: Feedback
    const id::Any
    soc::Real
    power_out::Real
    power_in::Real
end

mutable struct WindFeedback <: Feedback 
    const id::Any
    power::Real
    power_pot::Real
    v::Real
    r::Real
    q::Real
end

mutable struct DemandFeedback <: Feedback
    const id::Any
    power::Real
end

mutable struct SpotFeedback <: Feedback
    const id::Any
    buy::Real
    sell::Real
end

mutable struct SimpleMarketFeedback <: Feedback
    const id::Any
    power_buy::Real
    power_sell::Real
end