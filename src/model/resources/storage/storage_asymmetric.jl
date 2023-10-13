@doc raw"""
	storage_asymmetric!(EP::Model, inputs::Dict, setup::Dict)

Sets up variables and constraints specific to storage resources with asymmetric charge and discharge capacities. See ```storage()``` in ```storage.jl``` for description of constraints.
"""
function storage_asymmetric!(EP::Model, inputs::Dict, setup::Dict)
	# Set up additional variables, constraints, and expressions associated with storage resources with asymmetric charge & discharge capacity
	# (e.g. most chemical, thermal, and mechanical storage options with distinct charge & discharge components/processes)
	# STOR = 2 corresponds to storage with distinct power and energy capacity decisions and distinct charge and discharge power capacity decisions/ratings

	println("Storage Resources with Asmymetric Charge/Discharge Capacity Module")

	dfGen = inputs["dfGen"]
	Reserves = setup["Reserves"]
	CapacityReserveMargin = setup["CapacityReserveMargin"]

	T = inputs["T"]     # Number of time steps (hours)

	STOR_ASYMMETRIC = inputs["STOR_ASYMMETRIC"]

	### Constraints ###

	# Storage discharge and charge power (and reserve contribution) related constraints for symmetric storage resources:
	if Reserves == 1
		storage_asymmetric_reserves!(EP, inputs, setup)
	else
		if CapacityReserveMargin > 0
			# Maximum charging rate (including virtual charging to move energy held in reserve back to available storage) must be less than charge power rating
			@constraint(EP, [y in STOR_ASYMMETRIC, t in 1:T], EP[:vCHARGE][y,t] + EP[:vCAPRES_charge][y,t] <= EP[:eTotalCapCharge][y])
		else
			# Maximum charging rate (including virtual charging to move energy held in reserve back to available storage) must be less than charge power rating
			@constraint(EP, [y in STOR_ASYMMETRIC, t in 1:T], EP[:vCHARGE][y,t] <= EP[:eTotalCapCharge][y])
		end
	end

end

@doc raw"""
	storage_asymmetric_reserves!(EP::Model, inputs::Dict)

Sets up variables and constraints specific to storage resources with asymmetric charge and discharge capacities when reserves are modeled. See ```storage()``` in ```storage.jl``` for description of constraints.
"""
function storage_asymmetric_reserves!(EP::Model, inputs::Dict, setup::Dict)

	T = 1:inputs["T"]
	CapacityReserveMargin = setup["CapacityReserveMargin"] > 0

	STOR_ASYMMETRIC = inputs["STOR_ASYMMETRIC"]
	STOR_ASYM_REG = intersect(STOR_ASYMMETRIC, inputs["REG"]) # Set of asymmetric storage resources with REG reserves

    vCHARGE = EP[:vCHARGE]
    vREG_charge = EP[:vREG_charge]
    vCAPRES_charge = EP[:vCAPRES_charge]
    eTotalCapCharge = EP[:eTotalCapCharge]

    expr = @expression(EP, [y in STOR_ASYMMETRIC, t in T], 1 * vCHARGE[y, t]) # NOTE load-bearing "1 *"
    add_similar_to_expression!(expr[STOR_ASYM_REG, :], vREG_charge[STOR_ASYM_REG, :])
    if CapacityReserveMargin
        add_similar_to_expression!(expr[STOR_ASYMMETRIC, :], vCAPRES_charge[STOR_ASYMMETRIC, :])
    end
    @constraint(EP, [y in STOR_ASYMMETRIC, t in T], expr[y, t] <= eTotalCapCharge[y])
end
