@doc raw"""
    resources_with_maintenance(df::DataFrame)::Vector{Int}

    Get a vector of the R_ID's of all resources listed in a dataframe
    that have maintenance requirements. If there are none, return an empty vector.

    This method takes a specific dataframe because compound resources may have their
    data in multiple dataframes.
"""
function resources_with_fusion(df::DataFrame)::Vector{Int}
    if "FUS" in names(df)
        df[df.FUS.>0, :R_ID]
    else
        Vector{Int}[]
    end
end

@doc raw"""
    resources_with_maintenance(inputs::Dict)::Vector{Int}

    Get a vector of the R_ID's of all resources listed in a dataframe
    that have maintenance requirements. If there are none, return an empty vector.
"""
function resources_with_fusion(inputs::Dict)::Vector{Int}
    resources_with_fusion(inputs["dfTS"])
end

function fusion_average_net_electric_power_expression!(EP::Model, inputs::Dict)
    dfGen = inputs["dfGen"]

    G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)

    dfTS = inputs["dfTS"]
    by_rid(rid, sym) = by_rid_df(rid, sym, dfTS)

    FUS =  resources_with_fusion(inputs)

    #System-wide installed capacity is less than a specified maximum limit
    has_max_up = dfTS[dfTS.Max_Up .>= 0, :R_ID]
    has_max_up = intersect(has_max_up, FUS)

    active_frac = ones(G)
    avg_start_power = zeros(G)
    net_th_frac = ones(G)
    net_el_factor = zeros(G)

    active_frac[has_max_up] .= 1 .- by_rid(has_max_up,:Dwell_Time) ./ by_rid(has_max_up,:Max_Up)
    avg_start_power[has_max_up] .= by_rid(has_max_up,:Start_Energy) ./ by_rid(has_max_up,:Max_Up)
    net_th_frac[FUS] .= active_frac[FUS] .* (1 .- by_rid(FUS,:Recirc_Act)) .- by_rid(FUS,:Recirc_Pass) .- avg_start_power[FUS]
    net_el_factor[FUS] .= dfGen[FUS,:Eff_Down] .* net_th_frac[FUS]

    dfGen.Average_Net_Electric_Factor = net_el_factor

    @expression(EP, eCAvgNetElectric[y in FUS], EP[:vCCAP][y] * net_el_factor[y])
end

function fusion_systemwide_max_cap_constraint!(EP::Model, inputs::Dict)
    dfGen = inputs["dfGen"]
    dfTS = inputs["dfTS"]

    FUS =  resources_with_fusion(inputs)

    FIRST_ROW = 1
    col = :System_Max_Cap_MWe_net
    if string(col) in names(dfTS)
        max_cap = dfTS[FIRST_ROW, col]
        if max_cap >= 0
            @constraint(EP, cCSystemTot, sum(EP[:eCAvgNetElectric][FUS]) <= max_cap)
        end
    end
end

@doc raw"""
    fusion_constraints!(EP::Model, inputs::Dict)

Apply fusion-core-specific constraints to the model.

"""
function fusion_constraints!(EP::Model, inputs::Dict, setup::Dict)

    T = 1:inputs["T"]
    p = inputs["hours_per_subperiod"]
    ω = inputs["omega"]

    dfGen = inputs["dfGen"]
    dfTS = inputs["dfTS"]

    by_rid(rid, sym) = by_rid_df(rid, sym, dfTS)

    FUS = resources_with_fusion(inputs)
    vCP = EP[:vCP]
    vCCAP = EP[:vCCAP]
    vCSTART = EP[:vCSTART]
    vCCOMMIT = EP[:vCCOMMIT]
    core_cap_size(y) = by_rid(y, :Cap_Size)
    dwell_time(y) = by_rid(y, :Dwell_Time)
    max_starts(y) = by_rid(y, :Max_Starts)
    max_uptime(y) = by_rid(y, :Max_Up)
    recirc_passive(y) = by_rid(y, :Recirc_Pass)
    recirc_active(y) = by_rid(y, :Recirc_Act)
    start_energy(y) = by_rid(y, :Start_Energy)
    start_power(y) = by_rid(y, :Start_Power)

    eff_down(y) = dfGen[y, :Eff_Down]

    # Minimum and maximum core power output
    @constraints(EP, begin
        # Maximum thermal power generated by core y at hour y <= Max power of committed
        # core minus power lost from down time at startup
        [t in T, y in FUS], vCP[t,y] <= core_cap_size(y) * (vCCOMMIT[t,y] -
                                       dwell_time(y) * vCSTART[t,y])
    end)

    FINITE_STARTS = intersect(FUS, dfTS[dfTS.Max_Starts.>=0, :R_ID])

    #Limit on total core starts per year
    @constraint(EP, [y in FINITE_STARTS],
        sum(vCSTART[t,y] * ω[t] for t in T) <= max_starts(y) * vCCAP[y] / core_cap_size(y)
    )

    MAX_UPTIME = intersect(FUS, dfTS[dfTS.Max_Up.>=0, :R_ID])
    # TODO: throw error if Max_Up == 0 since it's confusing & illdefined

    # Core max uptime. If this parameter > 0,
    # the fusion core must be cycled at least every n hours.
    # Looks back over interior timesteps and ensures that a core cannot
    # be committed unless it has been started at some point in
    # the previous n timesteps
    @constraint(EP, [t in T, y in MAX_UPTIME],
            vCCOMMIT[t,y] <= sum(vCSTART[hoursbefore(p, t, 0:(max_uptime(y)-1)), y]))

    # Passive recirculating power, depending on built capacity
    @expression(EP, ePassiveRecircFus[t in T, y in FUS],
                vCCAP[y] * eff_down(y) * recirc_passive(y))

    # Active recirculating power, depending on committed capacity
    @expression(EP, eActiveRecircFus[t in T, y in FUS],
                core_cap_size(y) * eff_down(y) * recirc_active(y) *
        (vCCOMMIT[t,y] - vCSTART[t,y] * dwell_time(y))
    )
    # Startup energy, taken from the grid every time the core starts up
    @expression(EP, eStartEnergyFus[t in T, y in FUS],
                core_cap_size(y) * vCSTART[t,y] * eff_down(y) * start_energy(y))

    # Startup power, required margin on the grid when the core starts
    @expression(EP, eStartPowerFus[t in T, y in FUS],
                core_cap_size(y) * vCSTART[t,y] * eff_down(y) * start_power(y))

end

function total_fusion_power_balance_expressions!(EP::Model, inputs::Dict)
    T = 1:inputs["T"]     # Time steps
    Z = 1:inputs["Z"]     # Zones
    dfGen = inputs["dfGen"]
    FUS = resources_with_fusion(inputs)

    #Total recirculating power at each timestep
    @expression(EP, eTotalRecircFus[t in T, y in FUS],
                EP[:ePassiveRecircFus][t,y] + EP[:eActiveRecircFus][t,y] + EP[:eStartEnergyFus][t,y])

    # Total recirculating power from fusion in each zone
    gen_in_zone(z) = dfGen[dfGen.Zone .== z, :R_ID]

    FUS_IN_ZONE = [intersect(FUS, gen_in_zone(z)) for z in Z]
    @expression(EP, ePowerBalanceRecircFus[t in T, z in Z],
        -sum(eTotalRecircFus[t,y] for y in FUS_IN_ZONE[z]))

    EP[:ePowerBalance] += ePowerBalanceRecircFus
end
