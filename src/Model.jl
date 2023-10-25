#Function for calculating Jₘₐₓ at optimal (Jₘₐₓₒₚₜ) temperature from Nitrogen per leaf mass (Nₘ_f)
function Calc_Jₘₐₓ(Nₘ_f::S,a::T,b::T) where {S<:Real,T<:Float64}     
    return max(a*Nₘ_f+b,0.0)
end
#Function for calcualting the seasonal peak Jₘₐₓ from Nitrogen per leaf mass (Nₘ_f)
function Calc_Jₘₐₓ(Nₘ_f::S,a::T,b::T,b_Jₘₐₓ::T) where {S<:Real,T<:Float64}     
    return b_Jₘₐₓ*Calc_Jₘₐₓ(Nₘ_f,a,b)
end
#Function for calcualting Jₘₐₓ from Nitrogen per leaf mass (Nₘ_f)
function Calc_Jₘₐₓ(Nₘ_f::S,a::T,b::T,b_Jₘₐₓ::T,xₜ::T) where {S<:Real,T<:Float64}     
    return xₜ*Calc_Jₘₐₓ(Nₘ_f,a,b,b_Jₘₐₓ) 
end

#Function for calculating the quantum yield from Jmax
function Calc_α(xₜ::T,α_max::T) where {T<:Float64} 
    return xₜ*α_max
end

Calc_LAI(model::CCPHStruct) = model.treesize.Wf/model.treepar.LMA*model.treesize.N

#Calculate the irradiance incident on a leaf at irradiance I
Calc_Iᵢ(I::Float64,model::CCPHStruct) = I*model.treepar.k/(1-model.treepar.m)
Calc_Iᵢ(I::Float64,treepar::TreePar) = I*treepar.k/(1-treepar.m)

#calcualte total conductance
Calc_gₜ(gₛ::Real,model::CCPHStruct) = gₛ*model.treepar.r_gₛ
Calc_gₜ(gₛ::Real,treepar::TreePar) = gₛ*treepar.r_gₛ

#Calculate trait optimization objective function
function calc_gain(A::S,
    Jₘₐₓ::S,
    N_cost::T,
    E_cost::S) where {T<:Real,S<:Real}
    gain = A*E_cost-N_cost*Jₘₐₓ
    return gain
end

#Run the Coupled Canopy Photosynthesis and Hydraulics model
function CCPH_run(gₛ::S,Nₘ_f::S,growthlength::T,model::CCPHStruct) where {S<:Real,T<:Float64}    
    
    #Calculate Jₘₐₓ
    Jₘₐₓ = Calc_Jₘₐₓ(Nₘ_f,model.treepar.a_Jmax,model.treepar.b_Jmax,model.photopar.b_Jmax,model.treepar.Xₜ)    
    #Quantum yield
    model.photopar.α = Calc_α(model.treepar.Xₜ,model.treepar.α_max)
    #Irradiance incident on a leaf at canopy top
    Iᵢ = Calc_Iᵢ(model.env.I₀,model)
    #Calculate LAI 
    LAI = Calc_LAI(model)
    #calcualte total conductance
    gₜ = Calc_gₜ(gₛ,model)
    #calculate electron transport
    J = Calc_J(Iᵢ,Jₘₐₓ,model.photopar.α,model.photopar.θ)
    #Calcualte intercellular carbon dioxide concentration
    cᵢ = calc_opt_cᵢ(gₜ,J,model.photopar,model.env)    
    #Calculate leaf C assimilation
    A = calc_Assimilation(gₜ,cᵢ,model.env.P,model.env.Cₐ)
    #Calculate leaf transpiration
    E = calc_E(gₛ,model.env.VPD,model.env.P;cons=model.cons)
    #Calculate per tree carbon assimilation
    P = GPP(A,LAI,growthlength,model)
    #Calculate cost factor of hydraulic failure
    E_cost, Kₓₗ, ψ_c = Calc_K_cost(gₛ,model)
    #Carbon cost of nitrogen uptake and protein maintenance 
    N_cost = model.treepar.Nₛ
    #Calculate trait optimization objective function
    gain = calc_gain(A,Jₘₐₓ,N_cost,E_cost)
    αr = zero(eltype(P))
    modeloutput = CCPHOutput(P,αr,ψ_c,Kₓₗ,E_cost,gain,cᵢ,A,E)

    return modeloutput
end

function Trait_objective_fun(x::S,growthlength::T,model::CCPHStruct) where {S<:AbstractArray,T<:Float64} 
    gₛ,Nₘ_f = x    

    modeloutput = CCPH_run(gₛ,Nₘ_f,growthlength,model)

    return -modeloutput.Gain
end 

#Find optimal triats
function CCPHTraitmodel(growthlength::T,model::CCPHStruct;
    gₛ_guess::T=0.02,gₛ_lim_lo::T=0.001,gₛ_lim_hi::T=0.5,
    Nₘ_f_guess::T=0.012,Nₘ_f_lim_lo::T=0.001,Nₘ_f_lim_hi::T=0.05) where {T<:Float64}  
    
    x0 = [gₛ_guess, Nₘ_f_guess]    
  
    lower = [gₛ_lim_lo, Nₘ_f_lim_lo]   
    upper = [min(gₛ_lim_hi,0.5), Nₘ_f_lim_hi] 
    
    df = Optim.OnceDifferentiable(x->Trait_objective_fun(x,growthlength,model),x0;autodiff = :forward)   
    
    inner_optimizer = Optim.BFGS(linesearch = Optim.LineSearches.BackTracking())
    opt_trait = Optim.optimize(df, lower, upper, x0, Optim.Fminbox(inner_optimizer))
  
    Optim.converged(opt_trait)||error("No optimal traits could be found")    
    
    gₛ_opt = opt_trait.minimizer[1]
    Nₘ_f_opt = opt_trait.minimizer[2]    
    
    return gₛ_opt,Nₘ_f_opt
end

#initiate weather parameters from WeatherTS
function Init_weather_par!(i::Integer,model::CCPHStruct,weatherts::WeatherTS,photo_kinetic::PhotoKineticRates) 
    growthlength = weatherts.tot_annual_daylight[i] #Length of day light time during the gorwth period of a specific year (s)
    step_length = weatherts.daylight[i]/weatherts.tot_annual_daylight[i] #length of simulation step (proportion of a year)
    model.env.I₀ = weatherts.PAR[i] 
    model.env.VPD = weatherts.VPD[i]  
    model.env.θₛ = weatherts.θₛ[i]   
    model.env.Tₐ = weatherts.temp[i]
    PhotoPar!(model.photopar,photo_kinetic,weatherts.temp[i])
    model.treepar.Xₜ = weatherts.acclimation_fac[i] #Account for acclimation 
    model.treepar.rₘ = model.treepar.rₘ_ref*exp(log(model.treepar.Q₁₀_rₘ)/10*(model.env.Tₐ-model.treepar.T_rₘ_ref)) #Calcualte respiraiton rate at current air temperature
    
    return growthlength,step_length
end