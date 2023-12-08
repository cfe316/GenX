@doc raw"""
	write_capacity(path::AbstractString, inputs::Dict, setup::Dict, EP::Model))

Function for writing the diferent capacities for the different generation technologies (starting capacities or, existing capacities, retired capacities, and new-built capacities).
"""
function write_capacity(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
	# Capacity decisions
	dfGen = inputs["dfGen"]
	MultiStage = setup["MultiStage"]

    scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1
    energy_scale = scale_factor
    currency_per_energy_scale= scale_factor

	capdischarge = zeros(size(inputs["RESOURCES"]))
	for i in inputs["NEW_CAP"]
		if i in inputs["COMMIT"]
			capdischarge[i] = value(EP[:vCAP][i])*dfGen[!,:Cap_Size][i]
		else
			capdischarge[i] = value(EP[:vCAP][i])
		end
	end

	retcapdischarge = zeros(size(inputs["RESOURCES"]))
	for i in inputs["RET_CAP"]
		if i in inputs["COMMIT"]
			retcapdischarge[i] = first(value.(EP[:vRETCAP][i]))*dfGen[!,:Cap_Size][i]
		else
			retcapdischarge[i] = first(value.(EP[:vRETCAP][i]))
		end
	end

	capacity_constraint_dual = zeros(size(inputs["RESOURCES"]))
	if :Max_Cap_MW in propertynames(dfGen)
		for y in dfGen[dfGen.Max_Cap_MW.>0, :R_ID]
			capacity_constraint_dual[y] = -dual.(EP[:cMaxCap][y])
		end
	end

	capcharge = zeros(size(inputs["RESOURCES"]))
	retcapcharge = zeros(size(inputs["RESOURCES"]))
	existingcapcharge = zeros(size(inputs["RESOURCES"]))
	for i in inputs["STOR_ASYMMETRIC"]
		if i in inputs["NEW_CAP_CHARGE"]
			capcharge[i] = value(EP[:vCAPCHARGE][i])
		end
		if i in inputs["RET_CAP_CHARGE"]
			retcapcharge[i] = value(EP[:vRETCAPCHARGE][i])
		end
		existingcapcharge[i] = MultiStage == 1 ? value(EP[:vEXISTINGCAPCHARGE][i]) : dfGen[!,:Existing_Charge_Cap_MW][i]
	end

	capenergy = zeros(size(inputs["RESOURCES"]))
	retcapenergy = zeros(size(inputs["RESOURCES"]))
	existingcapenergy = zeros(size(inputs["RESOURCES"]))
	for i in inputs["STOR_ALL"]
		if i in inputs["NEW_CAP_ENERGY"]
			capenergy[i] = value(EP[:vCAPENERGY][i])
		end
		if i in inputs["RET_CAP_ENERGY"]
			retcapenergy[i] = value(EP[:vRETCAPENERGY][i])
		end
		existingcapenergy[i] = MultiStage == 1 ? value(EP[:vEXISTINGCAPENERGY][i]) :  dfGen[i,:Existing_Cap_MWh]
	end
	if !isempty(inputs["VRE_STOR"])
		for i in inputs["VS_STOR"]
			if i in inputs["NEW_CAP_STOR"]
				capenergy[i] = value(EP[:vCAPENERGY_VS][i])
			end
			if i in inputs["RET_CAP_STOR"]
				retcapenergy[i] = value(EP[:vRETCAPENERGY_VS][i])
			end
			existingcapenergy[i] = dfGen[i,:Existing_Cap_MWh] # multistage functionality doesn't exist yet for VRE-storage resources
		end
	end

	dfCap = DataFrame(
		Resource = inputs["RESOURCES"], Zone = dfGen[!,:Zone],
		StartCap = MultiStage == 1 ? value.(EP[:vEXISTINGCAP]) : dfGen[!,:Existing_Cap_MW],
		RetCap = retcapdischarge[:],
		NewCap = capdischarge[:],
		EndCap = value.(EP[:eTotalCap]),
		CapacityConstraintDual = capacity_constraint_dual[:],
		StartEnergyCap = existingcapenergy[:],
		RetEnergyCap = retcapenergy[:],
		NewEnergyCap = capenergy[:],
		EndEnergyCap = existingcapenergy[:] - retcapenergy[:] + capenergy[:],
		StartChargeCap = existingcapcharge[:],
		RetChargeCap = retcapcharge[:],
		NewChargeCap = capcharge[:],
		EndChargeCap = existingcapcharge[:] - retcapcharge[:] + capcharge[:]
	)

    dfCap.StartCap *= energy_scale
    dfCap.RetCap *= energy_scale
    dfCap.NewCap *= energy_scale
    dfCap.EndCap *= energy_scale
    dfCap.CapacityConstraintDual *= currency_per_energy_scale
    dfCap.StartEnergyCap *= energy_scale
    dfCap.RetEnergyCap *= energy_scale
    dfCap.NewEnergyCap *= energy_scale
    dfCap.EndEnergyCap *= energy_scale
    dfCap.StartChargeCap *= energy_scale
    dfCap.RetChargeCap *= energy_scale
    dfCap.NewChargeCap *= energy_scale
    dfCap.EndChargeCap *= energy_scale

	total = DataFrame(
			Resource = "Total", Zone = "n/a",
			StartCap = sum(dfCap[!,:StartCap]), RetCap = sum(dfCap[!,:RetCap]),
			NewCap = sum(dfCap[!,:NewCap]), EndCap = sum(dfCap[!,:EndCap]),
			CapacityConstraintDual = "n/a",
			StartEnergyCap = sum(dfCap[!,:StartEnergyCap]), RetEnergyCap = sum(dfCap[!,:RetEnergyCap]),
			NewEnergyCap = sum(dfCap[!,:NewEnergyCap]), EndEnergyCap = sum(dfCap[!,:EndEnergyCap]),
			StartChargeCap = sum(dfCap[!,:StartChargeCap]), RetChargeCap = sum(dfCap[!,:RetChargeCap]),
			NewChargeCap = sum(dfCap[!,:NewChargeCap]), EndChargeCap = sum(dfCap[!,:EndChargeCap])
		)

	dfCap = vcat(dfCap, total)
	CSV.write(joinpath(path, "capacity.csv"), dfCap)
	return dfCap
end
