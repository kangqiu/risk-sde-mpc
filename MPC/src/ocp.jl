
function build_mpc(mpc::AbstractMPC, nodes::NamedArray{Node}, verbose = false, solver = Ipopt.Optimizer)
    ocp = Model(solver)

    # set_attribute(ocp, "hsllib", HSL_jll.libhsl_path)
    # set_attribute(ocp, "linear_solver", "ma57") 
    # set_attribute(ocp, "ma57_automatic_scaling", "yes")

    if verbose == false
        set_silent(ocp)
    end
    
    define_variables!(ocp, nodes, mpc)
    define_parameters!(ocp, nodes, mpc)
    
    add_node_constraints!(ocp, nodes, mpc)

    add_energybalance!(ocp, collect(nodes), mpc)
    define_objective!(ocp, nodes, mpc)
    return ocp
end

function solve!(ocp::Model, almost_solved=true)
    optimize!(ocp)
    @assert is_solved_and_feasible(ocp; allow_local = true,  allow_almost = almost_solved)
    if termination_status == (ALMOST_LOCALLY_SOLVED)
        print("Warning: acceptable solution found")
    end
end



function add_node_constraints!(ocp::Model, nodes::NamedArray, mpc::AbstractMPC)

    for n in nodes.array
        model_constraints!(ocp, n, mpc)
        state_constraints!(ocp, n, mpc)
        input_constraints!(ocp, n, mpc)
    end
end

function define_objective!(ocp::Model, nodes::NamedArray, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩

    V = 0.0
    V += cost_balance(ocp, mpc)

    for n in nodes.array
        V += get_cost(ocp, n, mpc)/last(𝒩)
    end

   @objective(ocp, Min, V)
end


function define_objective!(ocp::Model, nodes::NamedArray, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω = 1:mpc.model.ω
    quantiles = mpc.model.quantiles

    if typeof(mpc.risk_measure) == Expectation
        Vωk = repeat(Any[0.0], length(Ω), last(𝒩))
        for ω in Ω
            for k in 𝒩
                Vωk[ω, k] += cost_balance(ocp, mpc, ω, k)
                for n in nodes.array
                    Vωk[ω, k] += get_cost(ocp, n, mpc, ω, k)
                end
            end
        end
        Vωk = Vωk/last(𝒩)
        Vω = repeat(Any[0.0], length(Ω))
        #all scenarios are equiprobable by construction
        for ω in Ω
            Vω[ω] = sum(Vωk[ω, :])
        end

        V = sum(1/last(Ω) * Vω)

    elseif typeof(mpc.risk_measure) == CVaR
        V = 0.0
        @variable(ocp, η >= 0)
        Vωk = repeat(Any[0.0], length(Ω), last(𝒩))
        Vcvar_ωk = repeat(Any[0.0], length(Ω), last(𝒩))
        for ω in Ω
            for k in 𝒩
                Vωk[ω, k] += cost_balance(ocp, mpc, ω, k)
                for n in nodes.array
                    if n isa SimpleMarket
                        Vcvar_ωk[ω, k] += get_cost(ocp, n, mpc, ω, k)
                    else 
                        Vωk[ω, k] += get_cost(ocp, n, mpc, ω, k)
                    end
                end
            end
        end

        @variable(ocp, Vω[Ω]) 
        @constraint(ocp, [ω ∈ Ω], Vω[ω] == sum(Vωk[ω, :])/last(𝒩))


        @variable(ocp, Vcvar_ω[Ω]) 
        @constraint(ocp, [ω ∈ Ω], Vcvar_ω[ω] == sum(Vcvar_ωk[ω, :])/last(𝒩))

        @variable(ocp, w[Ω])
        @constraint(ocp, [ω ∈ Ω], w[ω] >= 0)
        @constraint(ocp, [ω ∈ Ω], w[ω] >= Vcvar_ω[ω]-η)
        V += η + 1/last(Ω) * sum(w) * 1/(1-mpc.risk_measure.α/100)
        V += 1/last(Ω) * sum(Vω)
    end
   @objective(ocp, Min, V)
end
