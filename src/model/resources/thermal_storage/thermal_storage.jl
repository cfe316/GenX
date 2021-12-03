"""
GenX: An Configurable Capacity Expansion Model
Copyright (C) 2021,  Massachusetts Institute of Technology
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
A complete copy of the GNU General Public License v2 (GPLv2) is available
in LICENSE.txt.  Users uncompressing this from an archive may not have
received this license file.  If not, see <http://www.gnu.org/licenses/>.
"""

@doc raw"""
    thermal_storage(EP::Model, setup::Dict, inputs::Dict)

"""
function thermal_storage(EP::Model, setup::Dict, inputs::Dict)

	println("Thermal Storage Module")

	dfGen = inputs["dfGen"]

	G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones



	START_SUBPERIODS = inputs["START_SUBPERIODS"]
	INTERIOR_SUBPERIODS = inputs["INTERIOR_SUBPERIODS"]
	hours_per_subperiod = inputs["hours_per_subperiod"] #total number of hours per subperiod


	###load thermal storage inputs
	TS = inputs["TS"]
	dfTS = inputs["dfTS"]


	@variables(EP, begin
	#thermal core variables
	vCP[y in TS, t = 1:T] >= 0 		#thermal core power for resource y at timestep t
	vCCAP[y in TS] >= 0  			#thermal core capacity for resource y

	#thermal storage variables
	vTS[y in TS, t = 1:T] >= 0		#thermal storage state of charge for resource y at timestep t
	vTSCAP[y in TS] >= 0			#thermal storage energy capacity for resource y
	end)

	### THERMAL CORE CONSTRAINTS ###
	@constraint(EP, cCPMax[y in TS, t=1:T], vCP[y,t] <= vCCAP[y])	#core power output must be <= installed capacity
	@constraint(EP, cCCAPMax[y in TS], vCCAP[y] <= dfTS[dfTS.R_ID.==y,:Max_Cap_MW_th][])	#Total installed capacity is less than specified maximum limit

	#variable cost of core operation
	@expression(EP, eCVar_Core[y in TS, t=1:T], inputs["omega"][t]*dfTS[dfTS.R_ID.==y,:Var_OM_Cost_per_MWh_th][]*vCP[y,t]) 	#variable cost at timestep t for thermal core y
	@expression(EP, eTotalCVarCoreT[t=1:T], sum(eCVar_Core[y,t] for y in TS))	#variable cost from all thermal cores at timestep t
	@expression(EP, eTotalCVarCore, sum(eTotalCVarCoreT[t] for t in 1:T))	#total variable cost for all thermal cores
	EP[:eObj] += eTotalCVarCore


	#core investment costs
	@expression(EP, eCFixed_Core[y in TS], dfTS[dfTS.R_ID.==y,:Fixed_Cost_per_MW_th][]*vCCAP[y])	#fixed cost for thermal core y
	@expression(EP, eTotalCFixedCore, sum(eCFixed_Core[y] for y in TS))		#total fixed costs for all thermal cores
	EP[:eObj] += eTotalCFixedCore

	### THERMAL STORAGE CONSTRAINTS ###
	@constraint(EP, cTSMax[y in TS, t=1:T], vTS[y,t] <= vTSCAP[y])	#storage state of charge must be <= installed capacity

	@constraint(EP, cTSocBalInterior[t in INTERIOR_SUBPERIODS, y in TS], vTS[y,t] ==			#thermal state of charge balance for interior timesteps: (previous SOC) - (discharge to turbines) - (turbine startup energy use) + (core power output) - (self discharge)
		vTS[y,t-1]-(1/dfTS[dfTS.R_ID.==y,:Eff_Therm][]*EP[:vP][y,t])
		-(1/dfTS[dfTS.R_ID.==y,:Eff_Therm][]*dfGen[!,:Start_Fuel_MMBTU_per_MW][y]*dfGen[!,:Cap_Size][y]*EP[:vSTART][y,t])
		+(dfTS[dfTS.R_ID.==y,:Eff_Up][]*vCP[y,t])-(dfTS[dfTS.R_ID.==y,:Self_Disch][]*vTS[y,t-1]))

	if !(setup["OperationWrapping"] == 1 && setup["LongDurationStorage"] == 1) #thermal SOC balance for start timesteps if LDES is not enabled
		@constraint(EP, cTSoCBalStart[t in START_SUBPERIODS, y in TS], vTS[y,t] ==
			vTS[y,t+hours_per_subperiod-1]-(1/dfTS[dfTS.R_ID.==y,:Eff_Therm][]*EP[:vP][y,t])
			-(1/dfTS[dfTS.R_ID.==y,:Eff_Therm][]*dfGen[!,:Start_Fuel_MMBTU_per_MW][y]*dfGen[!,:Cap_Size][y]*EP[:vSTART][y,t])
			+(dfTS[dfTS.R_ID.==y,:Eff_Up][]*vCP[y,t])-(dfTS[dfTS.R_ID.==y,:Self_Disch][]*vTS[y,t+hours_per_subperiod-1]))
	end

	#constraints if LDES is active
	if setup["OperationWrapping"] == 1 && setup["LongDurationStorage"] == 1

		REP_PERIOD = inputs["REP_PERIOD"]  # Number of representative periods

		dfPeriodMap = inputs["Period_Map"] # Dataframe that maps modeled periods to representative periods
		NPeriods = size(inputs["Period_Map"])[1] # Number of modeled periods

		MODELED_PERIODS_INDEX = 1:NPeriods
		REP_PERIODS_INDEX = MODELED_PERIODS_INDEX[dfPeriodMap[!,:Rep_Period] .== MODELED_PERIODS_INDEX]

		@variable(EP, vTSOCw[y in TS, n in MODELED_PERIODS_INDEX] >= 0)

		# Build up in storage inventory over each representative period w
		# Build up inventory can be positive or negative
		@variable(EP, vdTSOC[y in TS, w=1:REP_PERIOD])
		# Note: tw_min = hours_per_subperiod*(w-1)+1; tw_max = hours_per_subperiod*w
		@constraint(EP, cThermSoCBalLongDurationStorageStart[w=1:REP_PERIOD, y in TS],
					    vTS[y,hours_per_subperiod*(w-1)+1] == (1-dfTS[dfTS.R_ID.==y,:Self_Disch][])*(vTS[y,hours_per_subperiod*w]-vdTSOC[y,w])
						-(1/dfTS[dfTS.R_ID.==y,:Eff_Therm][]*EP[:vP][y,hours_per_subperiod*(w-1)+1])
						-(1/dfTS[dfTS.R_ID.==y,:Eff_Therm][]*dfGen[!,:Start_Fuel_MMBTU_per_MW][y]*dfGen[!,:Cap_Size][y]*EP[:vSTART][y,hours_per_subperiod*(w-1)+1])
						+(dfTS[dfTS.R_ID.==y,:Eff_Up][]*vCP[y,hours_per_subperiod*(w-1)+1]))

		# Storage at beginning of period w = storage at beginning of period w-1 + storage built up in period w (after n representative periods)
		## Multiply storage build up term from prior period with corresponding weight
		@constraint(EP, cThermSoCBalLongDurationStorageInterior[y in TS, r in MODELED_PERIODS_INDEX[1:(end-1)]],
						vTSOCw[y,r+1] == vTSOCw[y,r] + vdTSOC[y,dfPeriodMap[!,:Rep_Period_Index][r]])

		## Last period is linked to first period
		@constraint(EP, cThermSoCBalLongDurationStorageEnd[y in TS, r in MODELED_PERIODS_INDEX[end]],
						vTSOCw[y,1] == vTSOCw[y,r] + vdTSOC[y,dfPeriodMap[!,:Rep_Period_Index][r]])

		# Storage at beginning of each modeled period cannot exceed installed energy capacity
		@constraint(EP, cThermSoCBalLongDurationStorageUpper[y in TS, r in MODELED_PERIODS_INDEX],
						vTSOCw[y,r] <= vTSCAP[y])

		# Initial storage level for representative periods must also adhere to sub-period storage inventory balance
		# Initial storage = Final storage - change in storage inventory across representative period
		@constraint(EP, cThermSoCBalLongDurationStorageSub[y in TS, r in REP_PERIODS_INDEX],
						vTSOCw[y,r] == vTS[y,hours_per_subperiod*dfPeriodMap[!,:Rep_Period_Index][r]] - vdTSOC[y,dfPeriodMap[!,:Rep_Period_Index][r]])
	end

	#thermal storage investment costs
	@expression(EP, eCFixed_TS[y in TS], dfTS[dfTS.R_ID.==y,:Fixed_Cost_per_MWh_th][]*vTSCAP[y]) #fixed costs for thermal storage y
	@expression(EP, eTotalCFixedTS, sum(eCFixed_TS[y] for y in TS))	#total fixed costs for all thermal storage
	EP[:eObj] += eTotalCFixedTS
	#thermal storage steam turbine constraints
	# JAS: TODO ask about this
	@constraint(EP, cTSGenLimit[y in TS, t=1:T], vCP[y,t] <= vTS[y,t])
	@constraint(EP, cTSStartUp[y in TS, t=1:T], vCP[y,t] <= vTS[y,t])

	### FUSION CONSTRAINTS ###
	FUS =  dfTS[dfTS.FUS.>=1,:R_ID]
	if !isempty(FUS) #use fusion constraints if thermal cores tagged 'FUS' are present

		#UC variables for the fusion core, analogous to standard UC
		@variables(EP, begin
			vFCOMMIT[y in FUS, t=1:T] >= 0
			vFSTART[y in FUS, t=1:T] >= 0
			vFSHUT[y in FUS, t=1:T] >= 0
		end)

		#upper bounds on core commitment/start/shut
		@constraints(EP, begin
			[y in FUS, t=1:T], vFCOMMIT[y,t] <= vCCAP[y]/dfTS[dfTS.R_ID.==y,:Cap_Size][]
			[y in FUS, t=1:T], vFSTART[y,t] <= vCCAP[y]/dfTS[dfTS.R_ID.==y,:Cap_Size][]
			[y in FUS, t=1:T], vFSHUT[y,t] <= vCCAP[y]/dfTS[dfTS.R_ID.==y,:Cap_Size][]
		end)

		# Commitment state constraint linking startup and shutdown decisions
		@constraints(EP, begin
			# For Start Hours, links first time step with last time step in subperiod
			[y in FUS, t in START_SUBPERIODS], vFCOMMIT[y,t] == vFCOMMIT[y,(t+hours_per_subperiod-1)] + vFSTART[y,t] - vFSHUT[y,t]
			# For all other hours, links commitment state in hour t with commitment state in prior hour + sum of start up and shut down in current hour
			[y in FUS, t in INTERIOR_SUBPERIODS], vFCOMMIT[y,t] == vFCOMMIT[y,t-1] + vFSTART[y,t] - vFSHUT[y,t]
		end)


		#minimum and maximum core power output
		@constraints(EP, begin
			# Minimum stable thermal power generated by core y at hour y >= Min power of committed core
			[y in TS, t=1:T], vCP[y,t] >= dfTS[dfTS.R_ID.==y,:Min_Power][]*dfTS[dfTS.R_ID.==y,:Cap_Size][]*vFCOMMIT[y,t]

			# Maximum thermal power generated by core y at hour y <= Max power of committed core minus power lost from down time at startup
			[y in TS, t=1:T], vCP[y,t] <= dfTS[dfTS.R_ID.==y,:Cap_Size][]*vFCOMMIT[y,t] - dfTS[dfTS.R_ID.==y,:Down_Time][]*dfTS[dfTS.R_ID.==y,:Cap_Size][]*vFSTART[y,t]
		end)

		#core max up time. If this parameter != -1, the fusion core must be cycled at least every n hours.
		for y in FUS
			Up_Time = dfTS[dfTS.R_ID.==y,:Up_Time][]

			if Up_Time >= 0
				F_Up_Time_HOURS = [] # Set of hours in the summation term of the maximum up time constraint for the first subperiod of each representative period
				for s in START_SUBPERIODS
					F_Up_Time_HOURS = union(F_Up_Time_HOURS, (s+1):(s+Up_Time-1))
				end
				@constraints(EP, begin
					# Looks back over interior timesteps, ensures that a core cannot be committed unless it has been started at some point in the previous n timesteps
					[t in setdiff(INTERIOR_SUBPERIODS,F_Up_Time_HOURS)], vFCOMMIT[y,t] <= sum(vFSTART[y,e] for e=(t+1-dfTS[dfTS.R_ID.==y,:Up_Time][]):t)

					#wraps up-time constraint around period ends
					[t in F_Up_Time_HOURS], vFCOMMIT[y,t] <= sum(vFSTART[y,e] for e=(t-((t%hours_per_subperiod)-1):t))+sum(vFSTART[y,e] for e=((t+hours_per_subperiod-(t%hours_per_subperiod))-(dfTS[dfTS.R_ID.==y,:Up_Time][]-(t%hours_per_subperiod))+1):(t+hours_per_subperiod-(t%hours_per_subperiod)))
					[t in START_SUBPERIODS], vFCOMMIT[y,t] <= vFSTART[y,t]+sum(vFSTART[y,e] for e=((t+hours_per_subperiod-1)-(dfTS[dfTS.R_ID.==y,:Up_Time][]-1)+1):(t+hours_per_subperiod-1))
				end)
			end
		end

		#passive and active recirculating power for each fusion generator
		#passive recirculating power, depending on built capacity
		@expression(EP, ePassiveRecircFus[y in FUS], vCCAP[y]*dfTS[dfTS.R_ID.==y,:Eff_Therm][]*dfTS[dfTS.R_ID.==y,:Recirc_Pass][])
		#active recirculating power, depending on committed capacity
		@expression(EP, eActiveRecircFus[y in FUS, t=1:T], dfTS[dfTS.R_ID.==y,:Cap_Size][]*vFCOMMIT[y,t]*dfTS[dfTS.R_ID.==y,:Eff_Therm][]*dfTS[dfTS.R_ID.==y,:Recirc_Act][])
		#startup energy, taken from the grid every time the core starts up
		@expression(EP, eStartEnergyFus[y in FUS, t=1:T], dfTS[dfTS.R_ID.==y,:Cap_Size][]*vFSTART[y,t]*dfTS[dfTS.R_ID.==y,:Eff_Therm][]*dfTS[dfTS.R_ID.==y,:Start_Energy][])

		#total recirculating power from fusion in each zone
		@expression(EP, ePowerBalanceRecircFus[t=1:T, z=1:Z],
			-sum((ePassiveRecircFus[y]+eActiveRecircFus[y,t]+eStartEnergyFus[y,t]) for y in intersect(FUS, dfTS[dfTS[!,:Zone].==z,:][!,:R_ID])))

		EP[:ePowerBalance] += ePowerBalanceRecircFus
	end

return EP
end
