module CartesianRuns

export Interval,
       CartesianRunIndices,
       shift,
       complement,
       expand

include("types.jl")
include("construction.jl")
include("common.jl")
include("expand.jl")
include("complement.jl")
include("intersect.jl")
include("setdiff.jl")
include("union.jl")

end
