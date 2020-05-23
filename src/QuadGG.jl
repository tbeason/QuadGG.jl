module QuadGG

using Printf


export adaptlob, adaptsim


mutable struct WarnFlag{b<:Bool}
    nowarn::b
end



include("adaptlob.jl")
include("adaptsim.jl")


end # module
