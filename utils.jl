function initialize_history()
    history = Dict(
    "feedback" => [], 
    "action" =>[],
    "demand" =>[], 
    "windpower_forecast" => [],
    "openloop" => [])
    return history
end

function update_history!(history::Dict, actions, feedback, demand, samples, ocp::Model)
    push!(history["action"], deepcopy(actions));
    push!(history["feedback"],deepcopy(feedback));
    push!(history["demand"], demand);
    push!(history["windpower_forecast"], samples);
    openloop_results = Dict()
    openloop_results["balance"] =  value.(ocp[:balance])
    openloop_results["demand"] = value.(ocp[:demand])
    openloop_results["price_buy"] = value.(ocp[:price_buy])
    openloop_results["price_sell"] = value.(ocp[:price_sell])
    openloop_results["u"] = value.(ocp[:u])
    openloop_results["x"] = value.(ocp[:x])
    openloop_results["σ"] = value.(ocp[:σ])
    push!(history["openloop"], openloop_results)
end