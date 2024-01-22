module CCPH
#Import packages
import Dates, Optim, SpecialFunctions, ForwardDiff

export Constants, EnvironmentStruct, PhotoKineticRates, EnvironmentFunStruct,
PhotoPar, PhotoPar!, TreeSize, TreePar, HydraulicsPar, CCPHStruct, CCPHOutput,
CCPHInstOutput, WeatherTS, CCPH_inst_run, CCPH_run!, CCPHOpt

#include external files
include("Auxiliary.jl")
include("ModelStructs.jl")
include("Upscaling.jl")
include("Photosynthesis.jl")
include("Hydraulic.jl")
include("Model.jl")
end