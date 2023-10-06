function resources_with_conventional_thermal_core(inputs::Dict)::Vector{Int}
    dfTS = inputs["dfTS"]
    if "FUSION" in names(dfTS)
        dfTS[dfTS.FUSION.==0,:R_ID]
    else
        inputs["THERM_STOR"]
    end
end

function resources_with_resistive_heating(inputs::Dict)::Vector{Int}
    dfTS = inputs["dfTS"]
    if "RH" in names(dfTS)
        dfTS[dfTS.RH.==1,:R_ID]
    else
        Vector{Int}[]
    end
end

function split_LDS_and_nonLDS(inputs::Dict)
    df = inputs["dfGen"]
    TS = inputs["THERM_STOR"]
    rep_periods = inputs["REP_PERIOD"]
    if rep_periods > 1
        LDS = df[df.LDS.==1, :R_ID]
        TS_and_LDS = intersect(TS, LDS)
        TS_and_nonLDS = intersect(TS, LDS)
    else
        TS_and_LDS = Int[]
        TS_and_nonLDS = TS
    end
    TS_and_LDS, TS_and_nonLDS
end

@doc raw"""
    thermal_storage(EP::Model, inputs::Dict, setup::Dict)

"""
function thermal_storage!(EP::Model, inputs::Dict, setup::Dict)

    @info "Thermal Storage Module"

    thermal_storage_base_variables!(EP, inputs)
    thermal_storage_core_commit!(EP, inputs, setup)

    thermal_storage_capacity_costs!(EP, inputs)
    thermal_storage_variable_costs!(EP, inputs)
    thermal_storage_quantity_constraints!(EP, inputs)
    thermal_storage_resistive_heating_power_balance!(EP, inputs)

    TS_and_LDS, TS_and_nonLDS = split_LDS_and_nonLDS(inputs)
    if !isempty(TS_and_LDS)
        thermal_storage_lds_constraints!(EP, inputs)
    end

    thermal_storage_capacity_ratio_constraints!(EP, inputs)
    thermal_storage_duration_constraints!(EP, inputs)

    ### CONVENTIONAL CORE CONSTRAINTS ###
    CONV = resources_with_conventional_thermal_core(inputs)

    if !isempty(CONV)
        conventional_thermal_core_effective_electric_power_expression!(EP, inputs)
        conventional_thermal_core_systemwide_max_cap_constraint!(EP, inputs)
        conventional_thermal_core_constraints!(EP, inputs, setup)
    end

    ### FUSION CONSTRAINTS ###
    FUSION = resources_with_fusion(inputs)

    if !isempty(FUSION)
        fusion_average_net_electric_power_expression!(EP, inputs)
        fusion_systemwide_max_cap_constraint!(EP, inputs)
        fusion_constraints!(EP, inputs, setup)
    end

    MAINTENANCE = resources_with_maintenance(inputs)
    if !isempty(MAINTENANCE)
        sanity_check_maintenance(MAINTENANCE, inputs)
        maintenance_constraints!(EP, inputs, setup)
    end

    if !isempty(intersect(MAINTENANCE, FUSION))
        maintenance_fusion_modification!(EP, inputs)
    end

    if setup["CapacityReserveMargin"] > 0
        thermal_storage_capacity_reserve_margin!(EP, inputs)
    end

    thermal_core_emissions!(EP, inputs)

    # must be run after maintenance
    total_fusion_power_balance_expressions!(EP, inputs)
    return nothing
end

function thermal_storage_component_capacity!(EP::Model,
                                             inputs::Dict,
                                             rids::Vector{Int},


function thermal_storage_base_variables!(EP::Model, inputs::Dict)
    T = 1:inputs["T"]
    TS = inputs["THERM_STOR"]
    RH = resources_with_resistive_heating(inputs)
    CAN_BUILD = inputs["NEW_CAP"]

    can_build_co = intersect(TS, CAN_BUILD)
    can_build_ts = intersect(TS, CAN_BUILD)
    can_build_rh = intersect(RH, CAN_BUILD)

    CAN_RETIRE = inputs["RET_CAP"]
    can_retire_ts = intersect(TS, CAN_RETIRE)
    can_retire_rh = intersect(RH, CAN_RETIRE)

    # each needs vNEWCAP_subcomp
    #            vRETCAP_subcomp
    #            eTotalCap_subcomp

    # if multistage
    #            vExisting_TS_subcomp

    @variables(EP, begin
        # Thermal core variables
        vCP[t in T, y in TS] >= 0      #thermal core power for resource y at timestep t
        vCCAP[y in TS] >= 0             #thermal core capacity for resource y

        # Thermal storage variables
        vTS[t in T, y in TS] >= 0      #thermal storage state of charge for resource y at timestep t
        vTSCAP[y in TS] >= 0            #thermal storage energy capacity for resource y

        # resistive heating variables
        vRH[t in T, y in RH] >= 0      #electrical energy from grid
        vRHCAP[y in RH] >= 0            #RH power capacity for resource
    end)
end

function thermal_storage_core_commit!(EP::Model, inputs::Dict, setup::Dict)
    dfTS = inputs["dfTS"]
    T = 1:inputs["T"]
    CONV = resources_with_conventional_thermal_core(inputs)
    FUSION =  resources_with_fusion(inputs)
    THERM_COMMIT = inputs["THERM_COMMIT"]
    p = inputs["hours_per_subperiod"]

    CONV_COMMIT = intersect(THERM_COMMIT, CONV)
    set = union(FUSION, CONV_COMMIT)

    by_rid(rid, sym) = by_rid_df(rid, sym, dfTS)

    if isempty(set)
        return
    end

    @variables(EP, begin
        vCCOMMIT[t in T, y in set] >= 0 #core commitment status
        vCSTART[t in T, y in set] >= 0 #core startup
        vCSHUT[t in T, y in set] >= 0 #core shutdown
    end)

    if setup["UCommit"] == 1 # Integer UC constraints
        for y in set
            set_integer.(vCCOMMIT[:, y])
            set_integer.(vCSTART[:, y])
            set_integer.(vCSHUT[:,y])
        end
    end

    # Upper bounds on core commitment/start/shut, and optional maintenance variables
    @constraints(EP, begin
        [t in T, y in set], vCCOMMIT[t,y] <= EP[:vCCAP][y] / by_rid(y,:Cap_Size)
        [t in T, y in set], vCSTART[t,y] <= EP[:vCCAP][y] / by_rid(y,:Cap_Size)
        [t in T, y in set], vCSHUT[t,y] <= EP[:vCCAP][y] / by_rid(y,:Cap_Size)
    end)

    # Minimum and maximum core power output
    @constraints(EP, begin
        # Minimum stable thermal power generated by core y at
        # hour y >= Min power of committed core
        [t in T, y in set], EP[:vCP][t,y] >= by_rid(y, :Min_Power) * by_rid(y, :Cap_Size) * vCCOMMIT[t,y]
        [t in T, y in set], EP[:vCP][t,y] <= by_rid(y, :Cap_Size) * vCCOMMIT[t,y]
    end)

    # Commitment state constraint linking startup and shutdown decisions (Constraint #4)
    @constraint(EP, [t in T, y in set],
        vCCOMMIT[t,y] == vCCOMMIT[hoursbefore(p,t,1), y] + vCSTART[t,y] - vCSHUT[t,y])
end

function thermal_storage_capacity_costs!(EP::Model, inputs::Dict)
    dfTS = inputs["dfTS"]
    TS = inputs["THERM__STOR"]
    RH = resources_with_resistive_heating(inputs)

    by_rid(rid, sym) = by_rid_df(rid, sym, dfTS)

    vCCAP = EP[:vCCAP]
    vTSCAP = EP[:vTSCAP]
    vRHCAP = EP[:vRHCAP]

    # Core investment costs
    # fixed cost for thermal core y
    @expression(EP, eCFixed_Core[y in TS], by_rid(y, :Fixed_Cost_per_MW_th) * vCCAP[y])
    # total fixed costs for all thermal cores
    @expression(EP, eTotalCFixedCore, sum(eCFixed_Core[y] for y in TS))
    EP[:eObj] += eTotalCFixedCore

    # Thermal storage investment costs
    # Fixed costs for thermal storage y
    @expression(EP, eCFixed_TS[y in TS], by_rid(y, :Fixed_Cost_per_MWh_th) * vTSCAP[y])
    # Total fixed costs for all thermal storage
    @expression(EP, eTotalCFixedTS, sum(eCFixed_TS[y] for y in TS))
    EP[:eObj] += eTotalCFixedTS

    # Resistive heating investment costs
    # Fixed costs for resource y
    @expression(EP, eCFixed_RH[y in RH], by_rid(y, :Fixed_Cost_per_MW_RH) * vRHCAP[y])
    # Total fixed costs for all resistive heating
    @expression(EP, eTotalCFixedRH, sum(eCFixed_RH[y] for y in RH))
    EP[:eObj] += eTotalCFixedRH
end

function thermal_storage_variable_costs!(EP::Model, inputs::Dict)
    dfTS = inputs["dfTS"]
    TS = inputs["THERM_STOR"]
    T = 1:inputs["T"]
    ω = inputs["omega"]

    by_rid(rid, sym) = by_rid_df(rid, sym, dfTS)
    vCP = EP[:vCP]

    # Variable cost of core operation
    # Variable cost at timestep t for thermal core y
    @expression(EP, eCVar_Core[t in T, y in TS], ω[t] * (by_rid(y, :Var_OM_Cost_per_MWh_th) + inputs["TS_C_Fuel_per_MWh"][y][t]) * vCP[t,y])
    # Variable cost from all thermal cores at timestep t)
    @expression(EP, eTotalCVarCoreT[t in T], sum(eCVar_Core[t,y] for y in TS))
    # Total variable cost for all thermal cores
    @expression(EP, eTotalCVarCore, sum(eTotalCVarCoreT[t] for t in T))
    EP[:eObj] += eTotalCVarCore
end

function thermal_storage_capacity_limits!(EP::Model, inputs::Dict)
    dfTS = inputs["dfTS"]
    by_rid(rid, sym) = by_rid_df(rid, sym, dfTS)

    vCCAP = EP[:vCCAP]
    # Total installed capacity is less than specified maximum limit
    those_with_max_cap = dfTS[dfTS.Max_Cap_MW_th.>=0, :R_ID]
    @constraint(EP, cCCAPMax[y in those_with_max_cap], vCCAP[y] <= by_rid(y, :Max_Cap_MW_th))
end


function thermal_storage_quantity_constraints!(EP::Model, inputs::Dict)
    dfGen = inputs["dfGen"]
    T = 1:inputs["T"]
    p = inputs["hours_per_subperiod"]
    TS = inputs["THERM_STOR"]
    RH = resources_with_resistive_heating(inputs)
    vP = EP[:vP]
    vCP = EP[:vCP]
    vCCAP = EP[:vCCAP]
    vTSCAP = EP[:vTSCAP]
    vRHCAP = EP[:vRHCAP]
    vTS = EP[:vTS]
    vRH = EP[:vRH]

    ### THERMAL CORE CONSTRAINTS ###
    # Core power output must be <= installed capacity, including hourly capacity factors
    @constraint(EP, cCPMax[t in T, y in TS], vCP[t,y] <= vCCAP[y]*inputs["pP_Max"][y,t])

    ### THERMAL STORAGE CONSTRAINTS ###
    # Storage state of charge must be <= installed capacity
    @constraint(EP, cTSMax[t in T, y in TS], vTS[t,y] <= vTSCAP[y])

    # thermal state of charge balance for interior timesteps:
    # (previous SOC) - (discharge to turbines) - (turbine startup energy use) + (core power output) - (self discharge)
    @expression(EP, eTSSoCBalRHS[t in T, y in TS],
        vTS[hoursbefore(p, t, 1), y]
        - (1 / dfGen[y, :Eff_Down] * vP[y,t])
        - (1 / dfGen[y, :Eff_Down] * dfGen[y, :Start_Fuel_MMBTU_per_MW] * dfGen[y,:Cap_Size] * EP[:vSTART][y,t])
        + (dfGen[y,:Eff_Up] * vCP[t,y])
        - (dfGen[y,:Self_Disch] * vTS[hoursbefore(p, t, 1), y]))

    for y in RH, t in T
        add_to_expression!(EP[:eTSSoCBalRHS][t,y], vRH[t,y])
    end

    @constraint(EP, cTSSoCBal[t in T, y in TS], vTS[t,y] == eTSSoCBalRHS[t,y])

    ### RESISTIVE HEATING ###
    # Capacity constraint for RH
    @constraint(EP, cRHMax[t in T, y in RH], vRH[t, y] <= vRHCAP[y])
end

function thermal_storage_resistive_heating_power_balance!(EP::Model, inputs::Dict)
    dfGen = inputs["dfGen"]
    T = 1:inputs["T"]
    Z = 1:inputs["Z"]
    RH = resources_with_resistive_heating(inputs)
    vRH = EP[:vRH]
    @expression(EP, ePowerBalanceRH[t in T, z in Z],
        - sum(vRH[t, y] for y in intersect(RH, dfGen[dfGen[!, :Zone].==z, :R_ID])))
    EP[:ePowerBalance] += ePowerBalanceRH
end


function thermal_storage_lds_constraints!(EP::Model, inputs::Dict)
    dfGen = inputs["dfGen"]
    p = inputs["hours_per_subperiod"]
    REP_PERIOD = inputs["REP_PERIOD"]
    dfPeriodMap = inputs["Period_Map"]

    TS_and_LDS, TS_and_nonLDS = split_LDS_and_nonLDS(inputs)
    nperiods = nrow(dfPeriodMap)

    MODELED_PERIODS_INDEX = 1:nperiods
    REP_PERIODS_INDEX = MODELED_PERIODS_INDEX[dfPeriodMap.Rep_Period .== MODELED_PERIODS_INDEX]

    vP = EP[:vP] # outflow
    vCP = EP[:vCP] # inflow
    vTS = EP[:vTS] # state of charge

    @variable(EP, vTSOCw[n in MODELED_PERIODS_INDEX, y in TS_and_LDS] >= 0)

    # Build up in storage inventory over each representative period w
    # Build up inventory can be positive or negative
    @variable(EP, vdTSOC[w in 1:REP_PERIOD, y in TS_and_LDS])
    # Note: tw_min = hours_per_subperiod*(w-1)+1; tw_max = hours_per_subperiod*w
    @constraint(EP, cThermSoCBalLongDurationStorageStart[w in 1:REP_PERIOD, y in TS_and_LDS], (
        vTS[hours_per_subperiod * (w - 1) + 1, y] ==
                   (1 - dfGen[y, :Self_Disch]) * (vTS[hours_per_subperiod * w, y] - vdTSOC[w, y])
                 - (1 / dfGen[y, :Eff_Down] * vP[y, hours_per_subperiod * (w - 1) + 1])
                 - (1 / dfGen[y, :Eff_Down] * dfGen[y,:Start_Fuel_MMBTU_per_MW] * dfGen[y,:Cap_Size] * EP[:vSTART][y,hours_per_subperiod * (w - 1) + 1])
                 + (dfGen[y, :Eff_Up] * vCP[y,hours_per_subperiod * (w - 1) + 1])
             ))

    # Storage at beginning of period w = storage at beginning of period w-1 + storage built up in period w (after n representative periods)
    ## Multiply storage build up term from prior period with corresponding weight
    @constraint(EP, cThermSoCBalLongDurationStorage[r in MODELED_PERIODS_INDEX, y in TS_and_LDS],
                    vTSOCw[mod1(r+1, nperiods), y] == vTSOCw[r, y] + vdTSOC[dfPeriodMap[r,:Rep_Period_Index], y])

    # Storage at beginning of each modeled period cannot exceed installed energy capacity
    @constraint(EP, cThermSoCBalLongDurationStorageUpper[r in MODELED_PERIODS_INDEX, y in TS_and_LDS],
                    vTSOCw[r, y] <= vTSCAP[y])

    # Initial storage level for representative periods must also adhere to sub-period storage inventory balance
    # Initial storage = Final storage - change in storage inventory across representative period
    @constraint(EP, cThermSoCBalLongDurationStorageSub[r in REP_PERIODS_INDEX, y in TS_and_LDS],
                    vTSOCw[r, y] == vTS[hours_per_subperiod*dfPeriodMap[r,:Rep_Period_Index], y]
                                    - vdTSOC[dfPeriodMap[r,:Rep_Period_Index], y])
end


function thermal_storage_capacity_ratio_constraints!(EP::Model, inputs::Dict)
    dfGen = inputs["dfGen"]
    dfTS = inputs["dfTS"]
    vCCAP = EP[:vCCAP]
    by_rid(rid, sym) = by_rid_df(rid, sym, dfTS)

    has_max_ratio = dfTS[dfTS.Max_Generator_Core_Power_Ratio.>=0, :R_ID]
    max_ratio(y) = by_rid(y, :Max_Generator_Core_Power_Ratio)
    @constraint(EP, cCPRatMax[y in has_max_ratio],
        vCCAP[y] * dfGen[y,:Eff_Down] * max_ratio(y) >= EP[:eTotalCap][y])

    has_min_ratio = dfTS[dfTS.Min_Generator_Core_Power_Ratio.>=0, :R_ID]
    min_ratio(y) = by_rid(y, :Min_Generator_Core_Power_Ratio)
    @constraint(EP, cCPRatMin[y in has_min_ratio],
        vCCAP[y] * dfGen[y,:Eff_Down] * min_ratio(y) <= EP[:eTotalCap][y])
end

function thermal_storage_duration_constraints!(EP::Model, inputs::Dict)
    dfGen = inputs["dfGen"]
    TS = inputs["THERM_STOR"]
    vCCAP = EP[:vCCAP]
    vTSCAP = EP[:vTSCAP]
    # Limits on storage duration
    MIN_DUR = intersect(TS, dfGen[dfGen.Min_Duration .>= 0, :R_ID])
    MAX_DUR = intersect(TS, dfGen[dfGen.Max_Duration .>= 0, :R_ID])
    @constraint(EP, cTSMinDur[y in MIN_DUR], vTSCAP[y] >= dfGen[y,:Min_Duration] * vCCAP[y])
    @constraint(EP, cTSMaxDur[y in MAX_DUR], vTSCAP[y] <= dfGen[y,:Max_Duration] * vCCAP[y])
end

function conventional_thermal_core_effective_electric_power_expression!(EP::Model, inputs::Dict)
    dfGen = inputs["dfGen"]

    # convert thermal capacities to electrical capacities
    CONV =  resources_with_conventional_thermal_core(inputs)
    @expression(EP, eCElectric[y in CONV], EP[:vCCAP][y] * dfGen[y, :Eff_Down])
end


function conventional_thermal_core_systemwide_max_cap_constraint!(EP::Model, inputs::Dict)
    dfTS = inputs["dfTS"]

    #System-wide installed capacity is less than a specified maximum limit
    FIRST_ROW = 1
    col = :Nonfus_System_Max_Cap_MWe
    if string(col) in names(dfTS)
        max_cap = dfTS[FIRST_ROW, col]
        if max_cap >= 0
            @constraint(EP, cNonfusSystemTot, sum(eCElectric[CONV]) <= max_cap)
        end
    end
end


# TODO make compatible with reserves
function conventional_thermal_core_constraints!(EP::Model, inputs::Dict, setup::Dict)
    CONV = resources_with_conventional_thermal_core(inputs)
    THERM_COMMIT = inputs["THERM_COMMIT"]
    THERM_NO_COMMIT = inputs["THERM_NO_COMMIT"]

    COMMIT = intersect(THERM_COMMIT, CONV)
    NON_COMMIT = intersect(THERM_NO_COMMIT, CONV)

    # constraints for generators not subject to UC
    if !isempty(NON_COMMIT)
        conventional_thermal_core_no_commit_constraints(EP, inputs)
    end

    # constraints for generatiors subject to UC
    if !isempty(COMMIT)
        conventional_thermal_core_commit_constraints(EP, inputs)
    end
end

function conventional_thermal_core_no_commit_constraints!(EP::Model, inputs::Dict)

    T = 1:inputs["T"]     # Number of time steps (hours)
    p = inputs["hours_per_subperiod"] #total number of hours per subperiod

    CONV = resources_with_conventional_thermal_core(inputs)
    THERM_NO_COMMIT = inputs["THERM_NO_COMMIT"]
    set = intersect(THERM_NO_COMMIT, CONV)

    vCP = EP[:vCP]
    vCCAP = EP[:vCCAP]

    dfTS = inputs["dfTS"]
    by_rid(rid, sym) = by_rid_df(rid, sym, dfTS)
    ramp_up_frac(y) = by_rid(y, :Ramp_Up_Frac)
    ramp_dn_frac(y) = by_rid(y, :Ramp_Dn_Frac)
    min_power(y) = by_rid(y, :Min_Power)

    # ramp up and ramp down rates
    @constraints(EP, begin
                     [t in T, y in set], vCP[t, y] - vCP[hoursbefore(p, t, 1), y] <= ramp_up_frac(y) * vCCAP[y]
                     [t in T, y in set], vCP[hoursbefore(p, t, 1), y] - vCP[t,y] <= ramp_dn_frac(y) * vCCAP[y]
    end)

    # minimum stable power
    @constraint(EP, [t in T, y in set], vCP[t,y] >= min_power(y) * vCCAP[y])
end

function conventional_thermal_core_commit_constraints!(EP::Model, inputs::Dict)

    T = 1:inputs["T"]     # Number of time steps (hours)
    p = inputs["hours_per_subperiod"] #total number of hours per subperiod
    ω = inputs["omega"]

    CONV = resources_with_conventional_thermal_core(inputs)
    THERM_COMMIT = inputs["THERM_COMMIT"]

    set = intersect(THERM_COMMIT, CONV)

    vCP = EP[:vCP]
    vCCAP = EP[:vCCAP]
    vCSTART = EP[:vCSTART]
    vCCOMMIT = EP[:vCCOMMIT]
    vCSHUT = EP[:vCSHUT]

    dfTS = inputs["dfTS"]
    by_rid(rid, sym) = by_rid_df(rid, sym, dfTS)
    cap_size(y) = by_rid(y, :Cap_Size)
    ramp_up_frac(y) = by_rid(y, :Ramp_Up_Frac)
    ramp_dn_frac(y) = by_rid(y, :Ramp_Dn_Frac)
    min_power(y) = by_rid(y, :Min_Power)

    up_time(y) = Int(floor(by_rid(y, :Up_Time)))
    down_time(y) = Int(floor(by_rid(y, :Down_Time)))

    ### Add startup costs ###
    @expression(EP, eCStartTS[t in T, y in set], (ω[t] * inputs["TS_C_Start"][y][t] * vCSTART[t, y]))
    @expression(EP, eTotalCStartTST[t in T], sum(eCStartTS[t,y] for y in set))
    @expression(EP, eTotalCStartTS, sum(eTotalCStartTST[t] for t in T))
    EP[:eObj] += eTotalCStartTS

    #ramp up
    @constraint(EP,[t in T, y in set],
                vCP[t,y]-vCP[hoursbefore(p, t, 1), y] <= ramp_up_frac(y)*cap_size(y)*(vCCOMMIT[t,y]-vCSTART[t,y])
                + min(1, max(min_power(y), ramp_up_frac(y)))*cap_size(y)*vCSTART[t,y]
                - min_power(y) * cap_size(y) * vCSHUT[t,y])

    #ramp down
    @constraint(EP,[t in T, y in set],
                vCP[hoursbefore(p, t, 1), y]-vCP[t,y] <= ramp_dn_frac(y)*cap_size(y)*(vCCOMMIT[t,y]-vCSTART[t,y])
                - min_power(y)*cap_size(y)*vCSTART[t,y]
                + min(1,max(min_power(y), ramp_dn_frac(y)))*cap_size(y)*vCSHUT[t,y])

    ### Minimum up and down times
    @constraint(EP, [t in T, y in set],
        vCCOMMIT[t,y] >= sum(vCSTART[hoursbefore(p, t, 0:(up_time(y) - 1)), y])
    )

    @constraint(EP, [t in T, y in set],
        vCCAP[y]/cap_size(y)-vCCOMMIT[t,y] >= sum(vCSHUT[hoursbefore(p, t, 0:(down_time(y) - 1)), y])
    )
end

function thermal_storage_capacity_reserve_margin!(EP::Model, inputs::Dict)
    dfGen = inputs["dfGen"]
    dfTS = inputs["dfTS"]
    T = 1:inputs["T"]
    reserves = 1:inputs["NCapacityReserveMargin"]
    capresfactor(res, y) = dfGen[y, Symbol("CapRes_$res")]

    TS = inputs["THERM_STOR"]
    FUSION = resources_with_fusion(inputs)
    CONV = resources_with_conventional_thermal_core(inputs)
    MAINTENANCE = get_maintenance(inputs)

    by_rid(rid, sym) = by_rid_df(rid, sym, dfTS)

    vP = EP[:vP]

    # @expression(EP, eCapResMarBalanceThermalStorageAdjustment[res in reserves, t in T],
    #             sum(capresfactor(res, y) * (vP[y,t] - EP[:eTotalCap][y]) for y in TS))

    # EP[:eCapResMarBalance] += eCapResMarBalanceThermalStorageAdjustment

    @expression(EP, eCapResMarBalanceFusionAdjustment[res in reserves, t in T],
                sum(capresfactor(res, y) * (- EP[:eStartPowerFus][t,y]
                                            - EP[:ePassiveRecircFus][t,y]
                                            - EP[:eActiveRecircFus][t,y]) for y in FUSION))

    EP[:eCapResMarBalance] += eCapResMarBalanceFusionAdjustment

    # remove plants from contributing while they are under maintenance
    FUSION_MAINT = intersect(FUSION, MAINTENANCE)
    if !isempty(FUSION_MAINT)
        avg_net_el_fus(y) = dfGen[y, :Average_Net_Electric_Factor] * by_rid(y, :Cap_Size)
        @expression(EP, eCapResMarBalanceFusionMaintAdj[res in reserves, t in T],
                    -sum(capresfactor(res, y) * EP[:vMDOWN][t, y] * avg_net_el_fus(y) for y in FUSION_MAINT))
        EP[:eCapResMarBalance] += eCapResMarBalanceFusionMaintAdj
    end

    CONV_MAINT = intersect(CONV, MAINTENANCE)
    if !isempty(CONV_MAINT)
        net_el_conv(y) = dfGen[y, :Eff_Down] * by_rid(y, :Cap_Size)
        @expression(EP, eCapResMarBalanceTSConvMaintAdj[res in reserves, t in T],
                    -sum(capresfactor(res, y) * EP[:vMDOWN][t, y] * net_el_conv(y) for y in CONV_MAINT))
        EP[:eCapResMarBalance] += eCapResMarBalanceTSConvMaintAdj
    end

end

function thermal_core_emissions!(EP::Model, inputs::Dict)
    dfTS = inputs["dfTS"]
    dfGen = inputs["dfGen"]

    TS = inputs["THERM_STOR"]   # R_IDs of resources with thermal storage
    G = 1:inputs["G"]
    T = 1:inputs["T"]
    Z = 1:inputs["Z"]

    CONV = resources_with_conventional_thermal_core(inputs)
    THERM_COMMIT = inputs["THERM_COMMIT"]
    by_rid(rid, sym) = by_rid_df(rid, sym, dfTS)

    @expression(EP, eEmissionsByPlantTS[y in G, t in T],
        if y ∉ TS
            0
        elseif y in intersect(THERM_COMMIT, CONV)
            by_rid(y, :CO2_per_MWh) * EP[:vCP][t, y] + by_rid(y, :CO2_per_Start) * EP[:vCSTART][t, y]
        else
            by_rid(y, :CO2_per_MWh) * EP[:vCP][t, y]
        end
    )

    @expression(EP, eEmissionsByZoneTS[z in Z, t in T], sum(eEmissionsByPlantTS[y,t] for y in intersect(TS, dfGen[(dfGen[!,:Zone].==z),:R_ID])))
        EP[:eEmissionsByPlant] += eEmissionsByPlantTS
        EP[:eEmissionsByZone] += eEmissionsByZoneTS
end

@doc raw"""
    maintenance_formulation_thermal_storage!(EP::Model, inputs::Dict, setup::Dict)

    Creates maintenance variables and constraints for thermal-commit plants.
"""
function maintenance_formulation_thermal_storage!(EP::Model, inputs::Dict, setup::Dict)

    @info "Maintenance Module for Thermal-Storage plants"

    ensure_maintenance_variable_records!(inputs)
    df = inputs["dfTS"]
    by_rid(rid, sym) = by_rid_df(rid, sym, df)

    MAINT = resources_with_maintenance(df)
    resource_component(y) = by_rid(y, :Resource)
    cap(y) = by_rid(y, :Cap_Size)
    maint_dur(y) = Int(floor(by_rid(y, :Maintenance_Duration)))
    maint_freq(y) = Int(floor(by_rid(y, :Maintenance_Cycle_Length_Years)))
    maint_begin_cadence(y) = Int(floor(by_rid(y, :Maintenance_Begin_Cadence)))

    integer_operational_unit_commitment = setup["UCommit"] == 1

    vcommit = :vCCOMMIT
    ecap = :vCCAP

    sanity_check_maintenance(MAINT, inputs)

    for y in MAINT
        maintenance_formulation!(EP,
                                inputs,
                                resource_component(y),
                                y,
                                maint_begin_cadence(y),
                                maint_dur(y),
                                maint_freq(y),
                                cap(y),
                                vcommit,
                                ecap,
                                integer_operational_unit_commitment)
    end
end
