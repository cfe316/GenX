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
    geothermal_ires(EP::Model, inputs::Dict, ires_inputs::Dict)

"""
function fusion(EP::Model, setup::Dict, inputs::Dict)

	dfTS = inputs["dfTS"]
	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones

	START_SUBPERIODS = inputs["START_SUBPERIODS"]
	INTERIOR_SUBPERIODS = inputs["INTERIOR_SUBPERIODS"]
	hours_per_subperiod = inputs["hours_per_subperiod"] #total number of hours per subperiod

	FUS =  dfTS[dfTS.FUS.>=1,:R_ID]


	#UC variables for the thermal core
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
		[y in FUS, t in START_SUBPERIODS], EP[:vCOMMIT][y,t] == EP[:vCOMMIT][y,(t+hours_per_subperiod-1)] + EP[:vSTART][y,t] - EP[:vSHUT][y,t]
		# For all other hours, links commitment state in hour t with commitment state in prior hour + sum of start up and shut down in current hour
		[y in FUS, t in INTERIOR_SUBPERIODS], EP[:vCOMMIT][y,t] == EP[:vCOMMIT][y,t-1] + EP[:vSTART][y,t] - EP[:vSHUT][y,t]
	end)

	## For Start Hours
	# Links last time step with first time step, ensuring position in hour 1 is within eligible ramp of final hour position
		# rampup constraints
	@constraint(EP,[y in FUS, t in START_SUBPERIODS],
		EP[:vP][y,t]-EP[:vP][y,(t+hours_per_subperiod-1)] <= dfGen[!,:Cap_Size][y]*EP[:vCOMMIT][y,t]
			- dfTS[dfTS.R_ID.==y,:Min_Power][]*dfGen[!,:Cap_Size][y]*EP[:vSHUT][y,t])

		# rampdown constraints
	@constraint(EP,[y in FUS, t in START_SUBPERIODS],
		EP[:vP][y,(t+hours_per_subperiod-1)]-EP[:vP][y,t] <= dfGen[!,:Cap_Size][y]*(EP[:vCOMMIT][y,t]-EP[:vSTART][y,t])
			- dfTS[dfTS.R_ID.==y,:Min_Power][]*dfGen[!,:Cap_Size][y]*EP[:vSTART][y,t]
			+ dfGen[!,:Cap_Size][y]*EP[:vSHUT][y,t])

	## For Interior Hours
		# rampup constraints
	@constraint(EP,[y in FUS, t in INTERIOR_SUBPERIODS],
		EP[:vP][y,t]-EP[:vP][y,t-1] <= dfGen[!,:Cap_Size][y]*EP[:vCOMMIT][y,t]
			- dfTS[dfTS.R_ID.==y,:Min_Power][]*dfGen[!,:Cap_Size][y]*EP[:vSHUT][y,t])

		# rampdown constraints
	@constraint(EP,[y in FUS, t in INTERIOR_SUBPERIODS],
		EP[:vP][y,t-1]-EP[:vP][y,t] <= dfGen[!,:Cap_Size][y]*(EP[:vCOMMIT][y,t]-EP[:vSTART][y,t])
			- dfTS[dfTS.R_ID.==y,:Min_Power][]*dfGen[!,:Cap_Size][y]*EP[:vSTART][y,t]
			+ dfGen[!,:Cap_Size][y]*EP[:vSHUT][y,t])


	@constraints(EP, begin
		# Minimum stable core power generated per technology "y" at hour "t" > Min power
		[y in TS, t=1:T], EP[:vP][y,t] >= dfTS[dfTS.R_ID.==y,:Min_Power][]*dfGen[!,:Cap_Size][y]*EP[:vCOMMIT][y,t]

		# Maximum core power generated per technology "y" at hour "t" < Max power
		[y in TS, t=1:T], EP[:vP][y,t] <= dfGen[!,:Cap_Size][y]*EP[:vCOMMIT][y,t]
	end)



### Expressions ###

## Power Balance Expressions ##


EP[:ePowerBalance] += ePowerBalanceIRES

## Objective Function Expressions ##



EP[:eObj] += eTotalCIRES


return EP
end
