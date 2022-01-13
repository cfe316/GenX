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
	corecappower = zeros(size(TSResources))
	for i in 1:TSG
		corecappower[i] = first(value.(EP[:vCCAP][dfTS[!,:R_ID][i]]))
	end

	corecapenergy = zeros(size(TSResources))
	for i in 1:TSG
		corecapenergy[i] = first(value.(EP[:vTSCAP][dfTS[!,:R_ID][i]]))
	end

	dfCoreCap = DataFrame(
		Resource = TSResources, Zone = dfTS[!,:Zone],
		CorePowerCap = corecappower[:],
		TSEnergyCap = corecapenergy[:]
	)

	if setup["ParameterScale"] ==1
		dfCoreCap.CorePowerCap = dfCoreCap.CorePowerCap * ModelScalingFactor
		dfCoreCap.TSEnergyCap = dfCoreCap.TSEnergyCap * ModelScalingFactor
	end
	CSV.write(joinpath(path,"TS_capacity.csv"), dfCoreCap)

	### CORE POWER TIME SERIES ###
	dfCorePwr = DataFrame(Resource = TSResources, Zone=dfTS[!,:Zone], AnnualSum = Array{Union{Missing,Float32}}(undef, TSG))
	P_therm = zeros(TSG,T)
	for i in 1:TSG
		if setup["ParameterScale"] ==1
			P_therm[i,:] = value.(EP[:vCP][dfTS[!,:R_ID][i],:]) * ModelScalingFactor
		else
			P_therm[i,:] = value.(EP[:vCP][dfTS[!,:R_ID][i],:])
		end
		dfCorePwr[!,:AnnualSum][i] = sum(inputs["omega"].* P_therm[i,:])
	end
	dfCorePwr = hcat(dfCorePwr, DataFrame(P_therm, :auto))
	auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("AnnualSum");[Symbol("t$t") for t in 1:T]]
	rename!(dfCorePwr,auxNew_Names)
	total = DataFrame(["Total" 0 sum(dfCorePwr[!,:AnnualSum]) fill(0.0, (1,T))], :auto)
	for t in 1:T
		total[:,t+3] .= sum(dfCorePwr[:,Symbol("t$t")][1:TSG])
	end
	rename!(total,auxNew_Names)
	dfCorePwr = vcat(dfCorePwr, total)
	CSV.write(joinpath(path,"TSCorePwr.csv"), dftranspose(dfCorePwr, false), writeheader=false)


	### THERMAL SOC TIME SERIES ###
	dfTSOC = DataFrame(Resource = TSResources, Zone=dfTS[!,:Zone], AnnualSum = Array{Union{Missing,Float32}}(undef, TSG))
	TSOC = zeros(TSG,T)
	for i in 1:TSG
		if setup["ParameterScale"] ==1
			TSOC[i,:] = value.(EP[:vTS][dfTS[!,:R_ID][i],:]) * ModelScalingFactor
		else
			TSOC[i,:] = value.(EP[:vTS][dfTS[!,:R_ID][i],:])
		end
		dfTSOC[!,:AnnualSum][i] = sum(inputs["omega"].* TSOC[i,:])
	end
	dfTSOC = hcat(dfTSOC, DataFrame(TSOC, :auto))
	auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("AnnualSum");[Symbol("t$t") for t in 1:T]]
	rename!(dfTSOC,auxNew_Names)
	total = DataFrame(["Total" 0 sum(dfTSOC[!,:AnnualSum]) fill(0.0, (1,T))], :auto)
	for t in 1:T
		total[:,t+3] .= sum(dfTSOC[:,Symbol("t$t")][1:TSG])
	end
	rename!(total,auxNew_Names)
	dfTSOC = vcat(dfTSOC, total)
	CSV.write(joinpath(path,"TS_SOC.csv"), dftranspose(dfTSOC, false), writeheader=false)

	### CORE STARTS TIMESERIES ###
	dfFStart = DataFrame(Resource = TSResources, Zone = dfTS[!,:Zone], Sum = Array{Union{Missing,Float32}}(undef, TSG))
	start = zeros(TSG,T)
	for i in TSG
		start[i,:] = value.(EP[:vFSTART])[dfTS[!,:R_ID][i],:]
		dfFStart[!,:Sum][i] = sum(start[i,:])
	end
	dfFStart = hcat(dfFStart, DataFrame(start, :auto))
	auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
	rename!(dfFStart,auxNew_Names)
	total = DataFrame(["Total" 0 sum(dfFStart[!,:Sum]) fill(0.0, (1,T))], :auto)
	for t in 1:T
		if v"1.3" <= VERSION < v"1.4"
			total[!,t+3] .= sum(dfFStart[:,Symbol("t$t")][1:TSG])
		elseif v"1.5" <= VERSION < v"1.7"
			total[:,t+3] .= sum(dfFStart[:,Symbol("t$t")][1:TSG])
		end
	end
	rename!(total,auxNew_Names)
	dfFStart = vcat(dfFStart, total)
	CSV.write(joinpath(path,"f_start.csv"), dftranspose(dfFStart, false), writeheader=false)

		### CORE SHUTS TIMESERIES ###
	dfFShut = DataFrame(Resource = TSResources, Zone = dfTS[!,:Zone], Sum = Array{Union{Missing,Float32}}(undef, TSG))
	shut = zeros(TSG,T)
	for i in TSG
		shut[i,:] = value.(EP[:vFSHUT])[dfTS[!,:R_ID][i],:]
		dfFShut[!,:Sum][i] = sum(shut[i,:])
	end
	dfFShut = hcat(dfFShut, DataFrame(shut, :auto))
	auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
	rename!(dfFShut,auxNew_Names)
	total = DataFrame(["Total" 0 sum(dfFShut[!,:Sum]) fill(0.0, (1,T))], :auto)
	for t in 1:T
		if v"1.3" <= VERSION < v"1.4"
			total[!,t+3] .= sum(dfFShut[:,Symbol("t$t")][1:TSG])
		elseif v"1.5" <= VERSION < v"1.7"
			total[:,t+3] .= sum(dfFShut[:,Symbol("t$t")][1:TSG])
		end
	end
	rename!(total,auxNew_Names)
	dfFShut = vcat(dfFShut, total)
	CSV.write(joinpath(path,"f_shut.csv"), dftranspose(dfFShut, false), writeheader=false)

    ### CORE STARTS TIMESERIES ###
    dfFCommit = DataFrame(Resource = TSResources, Zone = dfTS[!,:Zone], Sum = Array{Union{Missing,Float32}}(undef, TSG))
    commit = zeros(TSG,T)
    for i in TSG
        commit[i,:] = value.(EP[:vFCOMMIT])[dfTS[!,:R_ID][i],:]
        dfFCommit[!,:Sum][i] = sum(commit[i,:])
    end
    dfFCommit = hcat(dfFCommit, DataFrame(commit, :auto))
    auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
    rename!(dfFCommit,auxNew_Names)
    total = DataFrame(["Total" 0 sum(dfFCommit[!,:Sum]) fill(0.0, (1,T))], :auto)
    for t in 1:T
        if v"1.3" <= VERSION < v"1.4"
            total[!,t+3] .= sum(dfFCommit[:,Symbol("t$t")][1:TSG])
        elseif v"1.5" <= VERSION < v"1.7"
            total[:,t+3] .= sum(dfFCommit[:,Symbol("t$t")][1:TSG])
        end
    end
    rename!(total,auxNew_Names)
    dfFCommit = vcat(dfFCommit, total)
    CSV.write(joinpath(path,"f_commit.csv"), dftranspose(dfFCommit, false), writeheader=false)

	return dfCoreCap, dfCorePwr, dfTSOC, dfFStart, dfFShut, dfFCommit
end
