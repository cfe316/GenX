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
	write_costs(path::AbstractString, sep::AbstractString, inputs::Dict, setup::Dict, EP::Model)

Function for writing the costs pertaining to the objective function (fixed, variable O&M etc.).
"""
function write_fusion_cost_categories(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
	## Cost results
	dfGen = inputs["dfGen"]
	SEG = inputs["SEG"]  # Number of lines
	Z = inputs["Z"]     # Number of zones
	T = inputs["T"]     # Number of time steps (hours)

	msf = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1

	dfCostCategories = DataFrame(Costs = ["cFixVRE", "cFixStorage", "cFixOtherGen",  "cVar", "cNSE", "cSTART", "cUnmetRsv", "cNetworkExp"])
	cVar = (value(EP[:eTotalCVarOut])+ (!isempty(inputs["STOR_ALL"]) ? value(EP[:eTotalCVarIn]) : 0) + (!isempty(inputs["FLEX"]) ? value(EP[:eTotalCVarFlexIn]) : 0) - (!isempty(inputs["TS"]) ? sum(value(EP[:eCVar_out])[y,t] for y in inputs["TS"], t in 1:T) : 0)) * (msf^2)
	cFixVRE = (!isempty(inputs["VRE"]) ? sum(value(EP[:eCFix])[y] for y in inputs["VRE"]) : 0) * (msf^2)
	cFixStorage = ((!isempty(inputs["STOR_ALL"]) ? sum(value(EP[:eCFix])[y] for y in inputs["STOR_ALL"]) : 0)  + (!isempty(inputs["STOR_ALL"]) ? value(EP[:eTotalCFixEnergy]) : 0) + (!isempty(inputs["STOR_ASYMMETRIC"]) ? value(EP[:eTotalCFixCharge]) : 0)) * (msf^2)
	cFixOtherGen = (value(EP[:eTotalCFix]) - (!isempty(inputs["TS"]) ? sum(value(EP[:eCFix])[y] for y in inputs["TS"]) : 0)) * (msf^2) - cFixVRE - cFixStorage
	dfCost[!,Symbol("Total")] = [cFixVRE, cFixStorage, cFixOtherGen, cVar, value(EP[:eTotalCNSE]) * (msf^2), 0, 0, 0]

	if setup["UCommit"]>=1
		dfCost[!,2][6] = (value(EP[:eTotalCStart]) - (!isempty(inputs["TS"]) ? sum(value(EP[:eCStart])[y,t] for y in inputs["TS"], t in 1:T) : 0))  * (msf^2)
	end

	if setup["Reserves"]==1
		dfCost[!,2][7] = value(EP[:eTotalCRsvPen]) * (msf^2)
	end

	if setup["NetworkExpansion"] == 1 && Z > 1
		dfCost[!,2][8] = value(EP[:eTotalCNetworkExp]) * (msf^2)
	end
	CSV.write(joinpath(path,"fusion_cost_categories.csv"), dfCostCategories)
end
