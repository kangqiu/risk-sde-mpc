abstract type Node
end

abstract type Source <: Node
end

abstract type Sink <: Node
end

struct GasTurbine <: Source
    id::String
    ηₛ::Real                # nominal efficiency
    Fᵪ::Real                # copper losses
    LHV::Real               # lower heating value 
    Pₙₒₘ::Real              # nominal power output
    loadₘᵢₙ::Real           # minimum Load

    M_CH₄::Real             # molar mass natural gas
    M_CO₂::Real             # molar mass CO₂
end

function GasTurbine(id::String, Pₙₒₘ::Real, ηₛ::Real=0.38, Fᵪ::Real =0.43 , LHV::Real=50, loadₘᵢₙ::Real=0.0)

    return GasTurbine(
        id::String,
        ηₛ::Real, 
        Fᵪ::Real, 
        LHV::Real, 
        Pₙₒₘ::Real, 
        loadₘᵢₙ::Real,
        16.043,
        44.01)

end

struct Battery <: Node
    id::String
    rate::Real              # charge and discharge rate in MW
    capacity::Real          # energy capacity in MWs
    η::Real                 # efficiency
    socₘᵢₙ::Real
    socₘₐₓ::Real
end

function Battery(id::String, rate::Real, capacity::Real, η::Real)
    return(
        Battery(
            id::String, 
            rate::Real, 
            capacity::Real, 
            η::Real, 
            0.1, 
            0.9)
    )
end

struct WindTurbineSDE <: Source
    id::String
    Pₘₐₓ::Real              # nominal capacity
end


struct Demand <: Sink
    id::String
end

struct SimpleMarket <: Sink
    id::String
    capacity::Real #capacity in MW
end