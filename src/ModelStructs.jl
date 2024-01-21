#Structs used for Forest-Photosynthesis-Hydraulic model 

#Constants
mutable struct Constants{T<:Float64}
    M_H2O::T #Molar mass water (Kg mol⁻¹)
    M_C::T #Molar mass Carbon (Kg mol⁻¹
    ρ_H2O::T #Density water (Kg m⁻³)    
    g::T #gravitational acceleration (m s⁻²)
    r::T #ratio water conductance:co2 conductance (-)    
    ρ_vapor::T #Density of vapor (Kg m⁻³)
end
function Constants(;M_H2O::T=0.018,M_C::T=0.012,ρ_H2O::T=997.0,g::T=9.82,r::T=1.6,ρ_vapor::T=0.749) where {T<:Float64} 
    Constants(M_H2O,M_C,ρ_H2O,g,r,ρ_vapor)
end

#Struct containing  function to determine within-day environmental variables
#time t (time after sunrise); Time must be given in seconds
mutable struct EnvironmentFunStruct{F<:Function} 
    I₀::F #Above canopy irradiance (PAR mol m⁻² s⁻¹)
    Cₐ::F #Ambient carbon dioxide partial pressure (Pa)
    P::F #Atmospheric pressure (Pa)
    Tₐ::F #Ambient air tempreture (°C)
    VPD::F #Vapour-pressure deficit     
    ψₛ::F #Soil water potential (MPa)  
end

#Struct containing environmental variables
mutable struct EnvironmentStruct{T<:Float64} 
    I₀::T #Above canopy irradiance (PAR mol m⁻² s⁻¹)
    Cₐ::T #Ambient carbon dioxide partial pressure (Pa)
    P::T #Atmospheric pressure (Pa)
    Tₐ::T #Ambient air tempreture (°C)
    VPD::T #Vapour-pressure deficit     
    ψₛ::T #Soil water potential (MPa)  
end
#Standard vaules
EnvironmentStruct(;I₀::T=820.0*10^-6,Cₐ::T=400.0/10.0,P::T=1.0*10^5,
Tₐ::T=25.0,VPD::T=800.0,θₛ::T = 0.15) where {T<:Float64} =  EnvironmentStruct(I₀,Cₐ,P,Tₐ,VPD,θₛ)
function EnvironmentStruct!(t::Real,env::EnvironmentStruct,envfun::EnvironmentFunStruct)
    env.I₀ = envfun.I₀(t)
    env.Cₐ = envfun.Cₐ(t)
    env.P = envfun.P(t)
    env.Tₐ = envfun.Tₐ(t)
    env.VPD = envfun.VPD(t)
    env.ψₛ = envfun.ψₛ(t)
    return nothing
end

mutable struct ArrheniusKineticRate{T<:Float64}    
    K_ref::T #Kinetic rate at referance temperature 
    T_ref::T #Reference temperature (°C)
    Eₐ::T #Activation energy (J mol⁻¹)
end
#Calcualte Kinetic rate at temperature Temp (°C)
function (K::ArrheniusKineticRate)(T::Float64)
    T += 273.15 #Temperature in Kelvin
    T_ref = K.T_ref+273.15 #Reference temperature in Kelvin
    R = 8.3145 #Gas constant (J K⁻¹ mol⁻¹)
    return K.K_ref*exp(K.Eₐ*(1/T_ref-1/T)/R)
end

mutable struct EnzymeKineticRate{T<:Float64}
    Kₒₚₜ::T #Kinetic rate at optimal temperature 
    Tₒₚₜ::T #Optimal temperature (°C)
    E_a::T #Activation energy (J mol⁻¹)
    E_d::T #Deactivation energy (J mol⁻¹)       
end
#Calcualte Kinetic rate at temperature Temp (°C)
function (K::EnzymeKineticRate)(T::Float64)
    T += 273.15 #Temperature in Kelvin
    Tₒₚₜ = K.Tₒₚₜ+273.15 #Optimal temperature in Kelvin
    R = 8.3145 #Gas constant (J K⁻¹ mol⁻¹)
    return K.Kₒₚₜ*K.E_d*exp(K.E_a*(T-Tₒₚₜ)/(T*R*Tₒₚₜ))/(K.E_d-K.E_a*(1-exp(K.E_d*(T-Tₒₚₜ)/(T*R*Tₒₚₜ))))
end

mutable struct PhotoKineticRates{T<:ArrheniusKineticRate}
    K_C::T #Carboxylation Michaelis-Menten constant
    K_O::T #Oxygenation Michaelis-Menten constant
    Γ::T #The CO2 compensation point in the absence of mitochondrial (day) respiration (Pa)
    b_r::T #The rate of dark respiration to Vmax ratio (Pa)
    Jmax::EnzymeKineticRate #
    O::Float64 #Intercellular partial pressure of oxygen (Pa)
end
function PhotoKineticRates(;Temp_ref::T=25.0, #Reference temperature (°C)
    K_C_ref::T = 40.40, K_C_E_a::T = 59.36e3, #Reference value (Pa) and activation energy (J mol⁻¹) for carboxylation Michaelis-Menten constant 
    K_O_ref::T = 248.0e2, K_O_E_a::T = 35.94e3, #Reference value (Pa) and activation energy (J mol⁻¹) for oxygenation Michaelis-Menten constant
    Γ_ref::T = 4.17, Γ_E_a::T = 23.42e3, #Reference value (Pa) and activation energy (J mol⁻¹) for CO2 compensation point
    b_r_ref::T = 0.02, b_r_E_a::T = 7.88e3, #Reference value (-) and activation energy (J mol⁻¹) for dark respiration to Vmax ratio
    Jmaxₒₚₜ::T = 184.0e-6, Jmax_Tempₒₚₜ::T = 31.85, Jmax_E_a::T = 47.4e3, Jmax_E_d::T = 200.0e3,
    #Optimal value (mol m⁻² s⁻¹), optimal temperature (°C), activation (J mol⁻¹) and deactivaiton (J mol⁻¹) energy for Jamx
    O::T = 20.5e3) where {T<:Float64} 

    K_C = ArrheniusKineticRate(K_C_ref,Temp_ref,K_C_E_a)
    K_O = ArrheniusKineticRate(K_O_ref,Temp_ref,K_O_E_a)
    Γ = ArrheniusKineticRate(Γ_ref,Temp_ref,Γ_E_a)
    b_r = ArrheniusKineticRate(b_r_ref,Temp_ref,b_r_E_a)
    Jmax = EnzymeKineticRate(Jmaxₒₚₜ,Jmax_Tempₒₚₜ,Jmax_E_a,Jmax_E_d)    

    PhotoKineticRates(K_C,K_O,Γ,b_r,Jmax,O)
end

#Parameters values used in the Farquhar photosynthesis model
mutable struct PhotoPar{T<:Float64}
    K::T #Combined Michaelis-Menten constant (Pa)
    b_r::T #Ratio between rate of dark respiration and Vmax (Pa)
    Γ::T #The CO2 compensation point in the absence of mitochondrial (day) respiration (Pa)    
    α::T #quantum yield of electron tansport
    θ::T #Curvature of the light response curve
    b_Jmax::T #The ratio between Jmax at tempreatur Temp and Jmax at optimal tempreture (Jmaxₒₚₜ)    
end
#Reference values at temperature = 25°C
PhotoPar(;K::T=73.8,b_r::T=0.02,Γ::T=4.17,α::T=0.36,θ::T=0.7,b_Jmax::T=0.812) where {T<:Float64} = 
PhotoPar(K,b_r,Γ,α,θ,b_Jmax)
#Function for initiate PhotoPar object at a tempreture of Temp °C
function PhotoPar(kinetic::PhotoKineticRates,Temp::T;α::T=0.36,θ::T=0.7) where {T<:Float64}
    K = kinetic.K_C(Temp)*(1+kinetic.O/kinetic.K_O(Temp))
    b_r = kinetic.b_r(Temp)
    Γ = kinetic.Γ(Temp)
    b_Jmax = kinetic.Jmax(Temp)/kinetic.Jmax.Kₒₚₜ
    
    PhotoPar(K,b_r,Γ,α,θ,b_Jmax)
end
#Calc PhotoPar values at temperature Temp (⁰C)
function PhotoPar!(photo::PhotoPar,kinetic::PhotoKineticRates,Temp::Float64) 
    K = kinetic.K_C(Temp)*(1+kinetic.O/kinetic.K_O(Temp))
    b_r = kinetic.b_r(Temp)
    Γ = kinetic.Γ(Temp)
    b_Jmax = kinetic.Jmax(Temp)/kinetic.Jmax.Kₒₚₜ
    
    photo.K = K
    photo.b_r = b_r
    photo.Γ = Γ
    photo.b_Jmax = b_Jmax
        
    return nothing
end

#Struct containing variables which defines the size of the stand
mutable struct TreeSize{T<:Float64}       
    H::T #Tree height (to canopy top) (m)   
    LAI::T #Above ground leaf area index (m² m⁻²) 
    N::T #Stem density (# trees m⁻² ground area)    
end

#Struct containing tree and stand parameters
mutable struct TreePar{T<:Float64}   
    Nₛ::T #Carbon cost of nitrogen uptake and protein maintenance factor (-)    
    k::T #Light extinction coefficient (-)
    m::T #Average leaf transmittance (-)
    a_Jmax::T #Slope of the Nitrogen per leaf area (Nₐ)-Jmaxₒₚₜ line (mol m⁻² leaf s⁻¹ Nₐ⁻¹)
    b_Jmax::T #Intercept of the Nitrogen per leaf area (Nₐ)-Jmaxₒₚₜ line (mol m⁻² leaf s⁻¹)       
    r_gₛ::T #total leaf conductance (stomatal+mesophyll) to stomatal conductance (gₜ/gₛ)
    Xₜ::T #Factor [0,1] accounting for the delayed effect of temperature on gross primary production (-) 
    α_max::T #Seasonal maximum quantum yield (m² s mol)    
end
#Standard values
TreePar(;Nₛ::T=0.04,k::T=0.52,m::T=0.05,a_Jmax::T=0.033,b_Jmax::T=-1.1e-5,
r_gₛ::T = 0.42,Xₜ::T = 1.0,α_max::T = 0.36) where {T<:Float64}  = 
TreePar(Nₛ,k,m,a_Jmax,b_Jmax,r_gₛ,Xₜ,α_max)

#Parameteres used in the Hydraulics model
mutable struct HydraulicsPar{T<:Float64}
    ψ₅₀::T #water potential which causes 50% loss in leaf conductance (MPa)  
    b::T #Sensitivity of leaf conductance to water potential (-)
    i::T #Sensitivity to hydraulic failure
    Kₓₗ₀::T #Maximum xylem-leaf hydraulic conductance (mol m⁻² leaf s⁻¹ MPa⁻¹)     
end
#Standard values
HydraulicsPar(;ψ₅₀::T=-2.89,b::T=2.15,i::T=1.0,Kₓₗ₀::T=0.01) where {T<:Float64} = 
HydraulicsPar(ψ₅₀,b,i,Kₓₗ₀)

#Collection of structs used for the instantaneous Photosynthesis and Hydraulic model
mutable struct CCPHStruct
    cons::Constants
    env::EnvironmentStruct
    treepar::TreePar
    treesize::TreeSize   
    photopar::PhotoPar  
    hydPar::HydraulicsPar
end

#Struct collecting output from the instantaneous Coupled Canopy Photosynthesis and Hydraulic model
mutable struct CCPHInstOutput{T<:Real}    
    ψ_c::T #leaf water potential (MPa)
    Kₓₗ::T #root-to-canopy hydraulic conductance (mol m⁻² leaf s⁻¹ MPa⁻¹)
    K_cost::T #ratio between the root-canopy conductance and the maximal conductance (-), i.e., K_cost=Kₓₗ/Kₓₗ₀
    Gain::T #Instantaneous leaf performance measure (mol C m⁻² leaf area s⁻¹)
    cᵢ::T #intercellular carbon dioxide concentration (Pa)
    A::T #leaf C assimilation (mol C m⁻² leaf area s⁻¹)   
    E::T #leaf transpiration (mol H₂O m² leaf area s⁻¹)
    GPP::T #Above ground vegetation GPP (mol C m⁻² ground area s⁻¹) 
    Ec::T #Canopy transpiration  (mol H₂O m² ground area s⁻¹) 
end

mutable struct CCPHOutput{T<:Real} 
    Gain::T #Leaf performance measure (mol C m⁻² leaf area day⁻¹)
    GPP::T #Above ground vegetation GPP (g C m⁻² ground day⁻¹)
    Ec::T #Canopy transpiration (mm day⁻¹)
end

#Struct collecting times series results from growth simulation
mutable struct CCPHTS{T<:Float64}
    Wf::Array{T,1}
    Ww::Array{T,1}    
    H::Array{T,1}
    Hs::Array{T,1}
    As::Array{T,1}
    B::Array{T,1}
    N::Array{T,1}
    P::Array{T,1}
    αr::Array{T,1}
    ψ_c::Array{T,1}
    Kₓₗ::Array{T,1}
    K_cost::Array{T,1}
    Gain::Array{T,1}
    gₛ::Array{T,1}
    Nₘ_f::Array{T,1} 
    gₛ_crit::Array{T,1}   
end
CCPHTS() = CCPHTS(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],
Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])
function CCPHTS!(ccphts::CCPHTS,treesize::TreeSize)    
    push!(ccphts.Wf,treesize.Wf)
    push!(ccphts.Ww,treesize.Ww)
    push!(ccphts.H,treesize.H)
    push!(ccphts.Hs,treesize.Hs)
    push!(ccphts.As,treesize.As)
    push!(ccphts.B,treesize.B)
    push!(ccphts.N,treesize.N)
end
function CCPHTS!(ccphts::CCPHTS,modeloutput::CCPHOutput,gₛ::T,Nₘ_f::T,gₛ_crit::T) where {T<:Float64}
    push!(ccphts.P,modeloutput.P)
    push!(ccphts.αr,modeloutput.αr)
    push!(ccphts.ψ_c,modeloutput.ψ_c)
    push!(ccphts.Kₓₗ,modeloutput.Kₓₗ)
    push!(ccphts.K_cost,modeloutput.K_cost)
    push!(ccphts.Gain,modeloutput.Gain)
    push!(ccphts.gₛ,gₛ)
    push!(ccphts.Nₘ_f,Nₘ_f)
    push!(ccphts.gₛ_crit,gₛ_crit)
end

#Containter for weather related time series
mutable struct WeatherTS{T<:Float64}
    date::Array{Dates.DateTime,1} #Dates 
    daylight::Array{T,1} #Daylight lenght, for each timestep, given in seconds
    tot_annual_daylight::Array{T,1} #Total annual daylight lenght during growth period, given in seconds
    PAR::Array{T,1} #mol s⁻¹ m⁻²
    temp::Array{T,1} #Mean air temperature (⁰C)
    VPD::Array{T,1} #Pa
    θₛ::Array{T,1} #Volumetric soil water content 
    acclimation_fac::Array{T,1} #Factor [0,1] accounting for annual acclimation cycle, see Mäkelä 2004 & 2008  
end

function EnvironmentStruct(weatherts::WeatherTS,ind::Integer)
    I₀ = weatherts.PAR[ind] 
    VPD = weatherts.VPD[ind]  
    θₛ = weatherts.θₛ[ind]   
    Tₐ = weatherts.temp[ind]
    EnvironmentStruct(;I₀=I₀,Tₐ=Tₐ,VPD=VPD,θₛ=θₛ) 
end