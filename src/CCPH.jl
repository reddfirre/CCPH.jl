module CCPH
#Import packages
import Dates

export Constants, EnvironmentStruct, PhotoKineticRates, PhotoPar, PhotoPar!,
TreeSize, TreePar, HydraulicsPar, CCPHStruct, CCPHOutput, CCPHTS, WeatherTS, Simple_CCPH

#include external files
include("Auxiliary.jl")
include("ModelStructs.jl")
include("Photosynthesis.jl")
include("Hydraulic.jl")
include("Model.jl")

end