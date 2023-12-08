@doc raw"""
    load_minimum_capacity_requirement!(path::AbstractString, inputs::Dict, setup::Dict)

Read input parameters related to mimimum capacity requirement constraints (e.g. technology specific deployment mandates)
"""
function load_minimum_capacity_requirement!(path::AbstractString, inputs::Dict, setup::Dict)
    filename = "Minimum_capacity_requirement.csv"
    df = load_dataframe(joinpath(path, filename))
    inputs["NumberOfMinCapReqs"] = nrow(df)

    scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1

    inputs["MinCapReq"] = df[!, :Min_MW] / scale_factor # Convert to GW
    if hasproperty(df, :PriceCap)
        inputs["MinCapPriceCap"] = df[!, :PriceCap] / scale_factor
    end
    println(filename * " Successfully Read!")
end
