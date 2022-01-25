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
    thermal_storage(EP::Model, inputs::Dict)

"""
function thermal_storage(EP::Model, inputs::Dict)

	println("Thermal Storage Module")

	dfGen = inputs["dfGen"]

	G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones


	START_SUBPERIODS = inputs["START_SUBPERIODS"]
	INTERIOR_SUBPERIODS = inputs["INTERIOR_SUBPERIODS"]
	hours_per_subperiod = inputs["hours_per_subperiod"]


	# Load thermal storage inputs
	TS = inputs["TS"]
	dfTS = inputs["dfTS"]

	function by_rid(rid::Integer, sym::Symbol)
		return dfTS[dfTS.R_ID .== rid, sym][]
	end

	@variables(EP, begin
	# Thermal core variables
	vCP[y in TS, t = 1:T] >= 0 		#thermal core power for resource y at timestep t
	vCCAP[y in TS] >= 0  			#thermal core capacity for resource y

	# Thermal storage variables
	vTS[y in TS, t = 1:T] >= 0		#thermal storage state of charge for resource y at timestep t
	vTSCAP[y in TS] >= 0			#thermal storage energy capacity for resource y
	end)

	### THERMAL CORE CONSTRAINTS ###
	# Core power output must be <= installed capacity, including hourly capacity factors
	@constraint(EP, cCPMax[y in TS, t=1:T], vCP[y,t] <= vCCAP[y]*inputs["pP_Max"][y,t])
	# Total installed capacity is less than specified maximum limit
	those_with_max_cap = dfTS[dfTS.Max_Cap_MW_th.>0, :R_ID]
	@constraint(EP, cCCAPMax[y in those_with_max_cap], vCCAP[y] <= by_rid(y, :Max_Cap_MW_th))
	#System-wide installed capacity is less than a specified maximum limit
	FIRST_ROW = 1
	if dfTS[FIRST_ROW, :System_Max_Cap_MW_th] >= 0
		@constraint(EP, cCSystemTot, sum(vCCAP[y] for y in TS) == dfTS[FIRST_ROW, :System_Max_Cap_MW_th])
	end


	# Variable cost of core operation
	# Variable cost at timestep t for thermal core y
	@expression(EP, eCVar_Core[y in TS, t=1:T], inputs["omega"][t] * by_rid(y, :Var_OM_Cost_per_MWh_th) * vCP[y,t])
	# Variable cost from all thermal cores at timestep t)
	@expression(EP, eTotalCVarCoreT[t=1:T], sum(eCVar_Core[y,t] for y in TS))
	# Total variable cost for all thermal cores
	@expression(EP, eTotalCVarCore, sum(eTotalCVarCoreT[t] for t in 1:T))
	EP[:eObj] += eTotalCVarCore


	# Core investment costs
	# fixed cost for thermal core y
	@expression(EP, eCFixed_Core[y in TS], by_rid(y,:Fixed_Cost_per_MW_th) * vCCAP[y])
	# total fixed costs for all thermal cores
	@expression(EP, eTotalCFixedCore, sum(eCFixed_Core[y] for y in TS))
	EP[:eObj] += eTotalCFixedCore

	### THERMAL STORAGE CONSTRAINTS ###
	# Storage state of charge must be <= installed capacity
	@constraint(EP, cTSMax[y in TS, t=1:T], vTS[y,t] <= vTSCAP[y])

	# thermal state of charge balance for interior timesteps:
	# (previous SOC) - (discharge to turbines) - (turbine startup energy use) + (core power output) - (self discharge)
	@constraint(EP, cTSocBalInterior[t in INTERIOR_SUBPERIODS, y in TS], (
		vTS[y,t] == vTS[y,t-1]
		- (1 / dfGen[y, :Eff_Down] * EP[:vP][y,t])
		- (1 / dfGen[y, :Eff_Down] * dfGen[y, :Start_Fuel_MMBTU_per_MW] * dfGen[y,:Cap_Size] * EP[:vSTART][y,t])
		+ (dfGen[y,:Eff_Up] * vCP[y,t])
		- (dfGen[y,:Self_Disch] * vTS[y,t-1]))
	)

	# Thermal SOC balance constraints
	for y in TS
		#If LDS is not enabled
		if dfGen[y, :LDS] != 1
			@constraint(EP, cTSoCBalStart[t in START_SUBPERIODS, y],(
			 vTS[y,t] == vTS[y, t + hours_per_subperiod - 1]
				- (1 / dfGen[y, :Eff_Down] * EP[:vP][y,t])
				- (1 / dfGen[y, :Eff_Down] * dfGen[y, :Start_Fuel_MMBTU_per_MW] * dfGen[y, :Cap_Size] * EP[:vSTART][y,t])
				+ (dfGen[y, :Eff_Up] * vCP[y,t]) - (dfGen[y, :Self_Disch] * vTS[y,t + hours_per_subperiod - 1])
				))

		#If LDS is enabled
		else
			REP_PERIOD = inputs["REP_PERIOD"]  # Number of representative periods

			dfPeriodMap = inputs["Period_Map"] # Dataframe that maps modeled periods to representative periods
			NPeriods = size(inputs["Period_Map"])[1] # Number of modeled periods

			MODELED_PERIODS_INDEX = 1:NPeriods
			REP_PERIODS_INDEX = MODELED_PERIODS_INDEX[dfPeriodMap[!,:Rep_Period] .== MODELED_PERIODS_INDEX]

			@variable(EP, vTSOCw[y, n in MODELED_PERIODS_INDEX] >= 0)

			# Build up in storage inventory over each representative period w
			# Build up inventory can be positive or negative
			@variable(EP, vdTSOC[y, w=1:REP_PERIOD])
			# Note: tw_min = hours_per_subperiod*(w-1)+1; tw_max = hours_per_subperiod*w
			@constraint(EP, cThermSoCBalLongDurationStorageStart[w=1:REP_PERIOD, y], (
					vTS[y,hours_per_subperiod * (w - 1) + 1] ==
							   (1 - dfGen[y, :Self_Disch]) * (vTS[y, hours_per_subperiod * w] - vdTSOC[y,w])
							 - (1 / df_Gen[y, :Eff_Down] * EP[:vP][y, hours_per_subperiod * (w - 1) + 1])
							 - (1 / df_Gen[y, :Eff_Down] * dfGen[y,:Start_Fuel_MMBTU_per_MW] * dfGen[y,:Cap_Size] * EP[:vSTART][y,hours_per_subperiod * (w - 1) + 1])
						 + (dfGen[y, :Eff_Up] * vCP[y,hours_per_subperiod * (w - 1) + 1])
						 ))

			# Storage at beginning of period w = storage at beginning of period w-1 + storage built up in period w (after n representative periods)
			## Multiply storage build up term from prior period with corresponding weight
			@constraint(EP, cThermSoCBalLongDurationStorageInterior[y, r in MODELED_PERIODS_INDEX[1:(end-1)]],
							vTSOCw[y,r+1] == vTSOCw[y,r] + vdTSOC[y,dfPeriodMap[r,:Rep_Period_Index]])

			## Last period is linked to first period
			@constraint(EP, cThermSoCBalLongDurationStorageEnd[y, r in MODELED_PERIODS_INDEX[end]],
							vTSOCw[y,1] == vTSOCw[y,r] + vdTSOC[y,dfPeriodMap[r,:Rep_Period_Index]])

			# Storage at beginning of each modeled period cannot exceed installed energy capacity
			@constraint(EP, cThermSoCBalLongDurationStorageUpper[y, r in MODELED_PERIODS_INDEX],
							vTSOCw[y,r] <= vTSCAP[y])

			# Initial storage level for representative periods must also adhere to sub-period storage inventory balance
			# Initial storage = Final storage - change in storage inventory across representative period
			@constraint(EP, cThermSoCBalLongDurationStorageSub[y, r in REP_PERIODS_INDEX],
							vTSOCw[y,r] == vTS[y,hours_per_subperiod*dfPeriodMap[r,:Rep_Period_Index]] - vdTSOC[y,dfPeriodMap[r,:Rep_Period_Index]])

		end
	end

	# Thermal storage investment costs
	# Fixed costs for thermal storage y
	@expression(EP, eCFixed_TS[y in TS], by_rid(y,:Fixed_Cost_per_MWh_th) * vTSCAP[y])
	# Total fixed costs for all thermal storage
	@expression(EP, eTotalCFixedTS, sum(eCFixed_TS[y] for y in TS))
	EP[:eObj] += eTotalCFixedTS

	# Parameter Fixing Constraints
	# Fixed ratio of generator capacity to core net electric power
	@constraint(EP, cCPRat[y in dfTS[dfTS.Generator_Core_Power_Ratio.>0,:R_ID]],
				vCCAP[y] * dfGen[y,:Eff_Down] * (1 - by_rid(y,:Recirc_Pass) - by_rid(y,:Recirc_Act)) ==
				EP[:vCAP][y] * dfGen[y,:Cap_Size] / by_rid(y,:Generator_Core_Power_Ratio))
	# Fixed storage duration
	@constraint(EP, cTSDur[y in dfTS[dfTS.Duration.>0,:R_ID]], vTSCAP[y] == by_rid(y,:Duration) * vCCAP[y])

	### FUSION CONSTRAINTS ###
	FUS =  dfTS[dfTS.FUS.>=1,:R_ID]

	# TODO: define expressions & constraints for NONFUS
	NONFUS =  dfTS[dfTS.FUS.==0,:R_ID]

	# Use fusion constraints if thermal cores tagged 'FUS' are present
	if !isempty(FUS)

		# UC variables for the fusion core, analogous to standard UC
		@variables(EP, begin
			vFCOMMIT[y in FUS, t=1:T] >= 0
			vFSTART[y in FUS, t=1:T] >= 0
			vFSHUT[y in FUS, t=1:T] >= 0
		end)

		#Declare core integer/binary variables if Integer_Commit is set to 1
		for y in FUS
			if by_rid(y, :Integer_Commit) == 1
				set_integer.(vFCOMMIT[y,:])
				set_integer.(vFSTART[y,:])
				set_integer.(vFSHUT[y,:])
				set_integer.(vCCAP[y])
			end
		end

		# Upper bounds on core commitment/start/shut
		@constraints(EP, begin
			[y in FUS, t=1:T], vFCOMMIT[y,t] <= vCCAP[y] / by_rid(y,:Cap_Size)
			[y in FUS, t=1:T], vFSTART[y,t] <= vCCAP[y] / by_rid(y,:Cap_Size)
			[y in FUS, t=1:T], vFSHUT[y,t] <= vCCAP[y] / by_rid(y,:Cap_Size)
		end)

		# Commitment state constraint linking startup and shutdown decisions
		@constraints(EP, begin
			# For Start Hours, links first time step with last time step in subperiod
			[y in FUS, t in START_SUBPERIODS], vFCOMMIT[y,t] == vFCOMMIT[y,(t+hours_per_subperiod-1)] + vFSTART[y,t] - vFSHUT[y,t]
			# For all other hours, links commitment state in hour t with commitment state in
			# prior hour + sum of start up and shut down in current hour
			[y in FUS, t in INTERIOR_SUBPERIODS], vFCOMMIT[y,t] == vFCOMMIT[y,t-1] + vFSTART[y,t] - vFSHUT[y,t]
		end)


		# Minimum and maximum core power output
		@constraints(EP, begin
			# Minimum stable thermal power generated by core y at
			# hour y >= Min power of committed core
			[y in FUS, t=1:T], vCP[y,t] >= by_rid(y, :Min_Power) * by_rid(y, :Cap_Size) * vFCOMMIT[y,t]

			# Maximum thermal power generated by core y at hour y <= Max power of committed
			# core minus power lost from down time at startup
			[y in FUS, t=1:T], vCP[y,t] <= by_rid(y, :Cap_Size) * vFCOMMIT[y,t] -
			                              by_rid(y, :Dwell_Time) * by_rid(y, :Cap_Size) * vFSTART[y,t]
		end)


		for y in FUS
			#Limit on total core starts per year
			if by_rid(y, :Max_Starts) >= 0
				@constraints(EP, begin
					[t=1:T], sum(vFSTART[y,t]*inputs["omega"][t] for t in 1:T) <=
										by_rid(y, :Max_Starts) * vCCAP[y] / by_rid(y,:Cap_Size)
				end)
			end

			# Core max uptime. If this parameter > 0,
			# the fusion core must be cycled at least every n hours.
			Max_Up = by_rid(y, :Max_Up)

			# TODO: throw error if Max_Up == 0 since it's confusing & illdefined

			if Max_Up > 0
				# Set of hours in the summation term of the maximum up time constraint for
				# the first subperiod of each representative period
				F_Max_Up_HOURS = []
				for s in START_SUBPERIODS
					F_Max_Up_HOURS = union(F_Max_Up_HOURS, (s+1):(s+Max_Up-1))
				end
				@constraints(EP, begin
					# Looks back over interior timesteps and ensures that a core cannot
					# be committed unless it has been started at some point in
					# the previous n timesteps
					[t in setdiff(INTERIOR_SUBPERIODS,F_Max_Up_HOURS)], vFCOMMIT[y,t] <= sum(vFSTART[y,e] for e=(t+1-by_rid(y,:Max_Up)):t)

					# Wraps up-time constraint around period ends
					[t in F_Max_Up_HOURS], vFCOMMIT[y,t] <= sum(vFSTART[y,e] for e=(t-((t%hours_per_subperiod)-1):t)) +
														    sum(vFSTART[y,e] for e=((t+hours_per_subperiod-(t%hours_per_subperiod))-(by_rid(y,:Max_Up)-(t%hours_per_subperiod))+1):(t+hours_per_subperiod-(t%hours_per_subperiod)))
					[t in START_SUBPERIODS], vFCOMMIT[y,t] <= vFSTART[y,t] +
															  sum(vFSTART[y,e] for e=((t+hours_per_subperiod-1)-(by_rid(y,:Max_Up) - 1) + 1):(t+hours_per_subperiod-1))
				end)
			end
		end

		# Passive and active recirculating power for each fusion generator
		# Passive recirculating power, depending on built capacity
		@expression(EP, ePassiveRecircFus[y in FUS], vCCAP[y] * dfGen[y,:Eff_Down] * by_rid(y,:Recirc_Pass))
		# Active recirculating power, depending on committed capacity
		@expression(EP, eActiveRecircFus[y in FUS, t=1:T], by_rid(y,:Cap_Size) * vFCOMMIT[y,t] * dfGen[y,:Eff_Down] * by_rid(y,:Recirc_Act))
		# Startup energy, taken from the grid every time the core starts up
		@expression(EP, eStartEnergyFus[y in FUS, t=1:T], by_rid(y,:Cap_Size) * vFSTART[y,t] * dfGen[y,:Eff_Down] * by_rid(y,:Start_Energy))

		# Total recirculating power from fusion in each zone
		@expression(EP, ePowerBalanceRecircFus[t=1:T, z=1:Z],
			-sum((ePassiveRecircFus[y]+eActiveRecircFus[y,t]+eStartEnergyFus[y,t]) for y in intersect(FUS, dfTS[dfTS[!,:Zone].==z,:R_ID])))

		EP[:ePowerBalance] += ePowerBalanceRecircFus
	end

return EP
end
