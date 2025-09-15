struct SimulationParameters
    Δt::Int #in minutes
    start::DateTime
    stop::DateTime
    horizon::Int
end

function SimulationParameters(Δt::Int, start::DateTime, stop::DateTime)
    horizon = (stop - start)/Minute(Δt)
    return SimulationParameters(Δt, start, stop, horizon)
end

abstract type MPCParameters
end

struct MPCTimestruct <: MPCParameters
    Δt::Int  #timestep in minutes
    horizon::Int 
end