module CCPH
#Import packages
import Dates, Optim, DifferentialEquations, SpecialFunctions

export Constants, EnvironmentStruct, PhotoKineticRates, PhotoPar, PhotoPar!,
TreeSize, TreePar, HydraulicsPar, CCPHStruct, CCPHOutput, CCPHTS, WeatherTS, CCPH_run,
CCPHTraitmodel, CCPHStandGrowth!

#include external files
include("Auxiliary.jl")
include("ModelStructs.jl")
include("Photosynthesis.jl")
include("Hydraulic.jl")
include("Model.jl")
end