@doc raw"""
    load_maximum_capacity_requirement!(path::AbstractString, inputs::Dict, setup::Dict)

Read input parameters related to maximum capacity requirement constraints (e.g. technology specific deployment mandates)
"""
function load_maximum_capacity_requirement!(path::AbstractString, inputs::Dict, setup::Dict)
    filename = "Maximum_capacity_requirement.csv"
    df = load_dataframe(joinpath(path, filename))
    inputs["NumberOfMaxCapReqs"] = nrow(df)

    scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1

    inputs["MaxCapReq"] = df[!, :Max_MW] / scale_factor
    if hasproperty(df, :PriceCap)
        inputs["MaxCapPriceCap"] = df[!, :PriceCap] / scale_factor
    end
    println(filename * " Successfully Read!")
end
