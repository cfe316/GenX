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

function write_core_behaviors(symbol::Symbol, filename::AbstractString)
	df = DataFrame(Resource = TSResources, Zone = dfTS[!,:Zone], Sum = Array{Union{Missing,Float32}}(undef, TSG))
	event = zeros(TSG,T)
	for i in 1:TSG
		event[i,:] = value.(EP[symbol])[dfTS[i,:R_ID],:]
		df[i,:Sum] = sum(event[i,:])
	end
	df = hcat(df, DataFrame(event, :auto))
	auxNew_Names=[:Resource; :Zone; :Sum; [Symbol("t$t") for t in 1:T]]
	rename!(df,auxNew_Names)
	total = DataFrame(["Total" 0 sum(df[!,:Sum]) fill(0.0, (1,T))], :auto)
	for t in 1:T
		total[:,t+3] .= sum(df[:,Symbol("t$t")][1:TSG])
	end
	rename!(total,auxNew_Names)
	df = vcat(df, total)
	CSV.write(joinpath(path, filename), dftranspose(df, false), writeheader=false)

	return df
end

function write_scaled_values(symbol::Symbol, filename::AbstractString)
	df = DataFrame(Resource = TSResources, Zone=dfTS[!,:Zone], AnnualSum = Array{Union{Missing,Float32}}(undef, TSG))
	quantity = zeros(TSG,T)
	for i in 1:TSG
		quantity[i,:] = value.(EP[symbol][dfTS[i,:R_ID],:]) * msf
		df[i,:AnnualSum] = sum(inputs["omega"].* quantity[i,:])
	end
	df = hcat(df, DataFrame(quantity, :auto))
	auxNew_Names=[:Resource; :Zone; :AnnualSum; [Symbol("t$t") for t in 1:T]]
	rename!(df,auxNew_Names)
	total = DataFrame(["Total" 0 sum(df[!,:AnnualSum]) fill(0.0, (1,T))], :auto)
	for t in 1:T
		total[:,t+3] .= sum(df[:,Symbol("t$t")][1:TSG])
	end
	rename!(total,auxNew_Names)
	df = vcat(df, total)
	CSV.write(joinpath(path, filename), dftranspose(df, false), writeheader=false)

	return df
end

@doc raw"""
	write_capacity(path::AbstractString, inputs::Dict, setup::Dict, EP::Model))

Function for writing the diferent capacities for the different generation technologies (starting capacities or, existing capacities, retired capacities, and new-built capacities).
"""
function write_thermal_storage(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
	# Capacity decisions
	dfGen = inputs["dfGen"]
	dfTS = inputs["dfTS"]
	T = inputs["T"]


	TSResources = dfTS[!,:Resource]
	TSG = length(TSResources)
	corecappower = zeros(TSG)
	for i in 1:TSG
		corecappower[i] = first(value.(EP[:vCCAP][dfTS[i,:R_ID]]))
	end

	corecapenergy = zeros(TSG)
	for i in 1:TSG
		corecapenergy[i] = first(value.(EP[:vTSCAP][dfTS[i,:R_ID]]))
	end

	dfCoreCap = DataFrame(
		Resource = TSResources, Zone = dfTS[!,:Zone],
		CorePowerCap = corecappower[:],
		TSEnergyCap = corecapenergy[:]
	)

	# set a single scalar to avoid future branching
	msf = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1

	dfCoreCap.CorePowerCap = dfCoreCap.CorePowerCap * msf
	dfCoreCap.TSEnergyCap = dfCoreCap.TSEnergyCap * msf
	CSV.write(joinpath(path,"TS_capacity.csv"), dfCoreCap)

	### CORE POWER TIME SERIES ###
	dfCorePwr = write_scaled_values(:vCP, "TSCorePwr.csv")

	### THERMAL SOC TIME SERIES ###
	dfTSOC = write_scaled_values(:vTS, "TS_SOC.csv")

	### CORE STARTS, SHUTS, and COMMITS TIMESERIES ###
	dfFStart = write_core_behaviors(:vFSTART, "f_start.csv")
	dfFShut = write_core_behaviors(:vFSHUT, "f_shut.csv")
	dfFCommit = write_core_behaviors(:vFCOMMIT, "f_commit.csv")

	return dfCoreCap, dfCorePwr, dfTSOC, dfFStart, dfFShut, dfFCommit
end
