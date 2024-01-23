#Function for calculating Jₘₐₓ at optimal (Jₘₐₓₒₚₜ) temperature from Nitrogen per leaf mass (Nₘ_f)
function calc_Jₘₐₓ(Nₘ_f::S,a::T,b::T) where {S<:Real,T<:Float64}     
    return max(a*Nₘ_f+b,0.0)
end
#Function for calcualting the seasonal peak Jₘₐₓ from Nitrogen per leaf mass (Nₘ_f)
function calc_Jₘₐₓ(Nₘ_f::S,a::T,b::T,b_Jₘₐₓ::T) where {S<:Real,T<:Float64}     
    return b_Jₘₐₓ*calc_Jₘₐₓ(Nₘ_f,a,b)
end
#Function for calcualting Jₘₐₓ from Nitrogen per leaf mass (Nₘ_f)
function calc_Jₘₐₓ(Nₘ_f::S,a::T,b::T,b_Jₘₐₓ::T,xₜ::T) where {S<:Real,T<:Float64}     
    return xₜ*calc_Jₘₐₓ(Nₘ_f,a,b,b_Jₘₐₓ) 
end

#Function for calculating the quantum yield from Jmax
function calc_α(xₜ::T,α_max::T) where {T<:Float64} 
    return xₜ*α_max
end

#Calculate the irradiance incident on a leaf at irradiance I
calc_Iᵢ(I::Float64,model::CCPHStruct) = I*model.treepar.k/(1-model.treepar.m)
calc_Iᵢ(I::Float64,treepar::TreePar) = I*treepar.k/(1-treepar.m)

#calcualte total conductance
calc_gₜ(gₛ::Real,model::CCPHStruct) = gₛ*model.treepar.r_gₛ
calc_gₜ(gₛ::Real,treepar::TreePar) = gₛ*treepar.r_gₛ

#Calculate instantaneous reward
function calc_gain(A::S,
    Jₘₐₓ::S,
    N_cost::T,
    E_cost::S) where {T<:Real,S<:Real}
    gain = A*E_cost-N_cost*Jₘₐₓ
    return gain
end

#Run the instantaneous Coupled Canopy Photosynthesis and Hydraulics model
function CCPH_inst_run(gₛ::S,Nₘ_f::S,model::CCPHStruct) where {S<:Real}    
    
    #Calculate Jₘₐₓ
    Jₘₐₓ = calc_Jₘₐₓ(Nₘ_f,model.treepar.a_Jmax,model.treepar.b_Jmax,model.photopar.b_Jmax,model.treepar.Xₜ)    
    #Quantum yield
    model.photopar.α = calc_α(model.treepar.Xₜ,model.treepar.α_max)
    #Irradiance incident on a leaf at canopy top
    Iᵢ = calc_Iᵢ(model.env.I₀,model)    
    #calcualte total conductance
    gₜ = calc_gₜ(gₛ,model)
    #calculate electron transport
    J = calc_J(Iᵢ,Jₘₐₓ,model.photopar.α,model.photopar.θ)
    #Calcualte intercellular carbon dioxide concentration
    cᵢ = calc_opt_cᵢ(gₜ,J,model.photopar,model.env)    
    #Calculate leaf C assimilation
    A = calc_Assimilation(gₜ,cᵢ,model.env.P,model.env.Cₐ)
    #Calculate leaf transpiration
    E = calc_E(gₛ,model.env.VPD,model.env.P;cons=model.cons)    
    #Calculate cost factor of hydraulic failure
    E_cost, Kₓₗ, ψ_c = calc_K_cost(gₛ,model)
    #Carbon cost of nitrogen uptake and protein maintenance 
    N_cost = model.treepar.Nₛ
    #Calculate trait optimization objective function
    gain = calc_gain(A,Jₘₐₓ,N_cost,E_cost)
    #Calculate canopy transpiration 
    Ec = calc_Ec(E,model)
    #Calculate above ground vegetation GPP
    GPP = calc_GPP(A,model)
    
    modelinstoutput = CCPHInstOutput(ψ_c,Kₓₗ,E_cost,gain,cᵢ,A,E,GPP,Ec)

    return modelinstoutput
end

#Run the Coupled Canopy Photosynthesis and Hydraulics model to get daily output values
function CCPH_run!(gₛ₁::S,gₛ₂::S,Nₘ_f::S,daylength::T,photo_kinetic::PhotoKineticRates,envfun::EnvironmentFunStruct,model::CCPHStruct) where {S<:Real,T<:Real}
    #SDM-2 is used to approximate the time integration from sunrise to sunset (Wang et al. 2014) 
    
    t₁ = daylength*asin(2/π)/(2*π)
    t₂ = daylength/4+t₁
    Δt₁ = t₁*2
    Δt₂ = t₂*2

    Init_weather_par!(t₁,model,photo_kinetic,envfun)
    modelinstoutput₁ = CCPH_inst_run(gₛ₁,Nₘ_f,model)
    
    Init_weather_par!(t₂,model,photo_kinetic,envfun)
    modelinstoutput₂ = CCPH_inst_run(gₛ₂,Nₘ_f,model)

    #Calculate above ground vegetation GPP (mol C m⁻² ground area day⁻¹)
    GPP = 2*(modelinstoutput₁.GPP*Δt₁+modelinstoutput₂.GPP*Δt₂)
    GPP *= model.cons.M_C*1000 #(g C m⁻² ground area day⁻¹)

    #Calculate canopy transpiration (mol H₂O m⁻² ground day⁻¹)
    Ec = 2*(modelinstoutput₁.Ec*Δt₁+modelinstoutput₂.Ec*Δt₂)
    Ec *= 1000*model.cons.M_H2O/model.cons.ρ_H2O #(mm day⁻¹)   

    #Calculate Leaf performance measure (mol C m⁻² leaf area day⁻¹)
    Gain = 2*(modelinstoutput₁.Gain*Δt₁+modelinstoutput₂.Gain*Δt₂)

    modeloutput = CCPHOutput(Gain,GPP,Ec)

    return modeloutput
end   

function Objective_fun(x::S,daylength::T,photo_kinetic::PhotoKineticRates,envfun::EnvironmentFunStruct,model::CCPHStruct) where {S<:AbstractArray,T<:Real} 
    gₛ₁,gₛ₂,Nₘ_f = x 
    
    modeloutput = CCPH_run!(gₛ₁,gₛ₂,Nₘ_f,daylength,photo_kinetic,envfun,model)

    return -modeloutput.Gain
end 

#Find optimal gₛ and Nₘ_f
function CCPHOpt(daylength::Real,photo_kinetic::PhotoKineticRates,envfun::EnvironmentFunStruct,model::CCPHStruct;
    gₛ₁_guess::Real=0.02,gₛ₁_lim_lo::Real=0.001,gₛ₁_lim_hi::Real=0.5,
    gₛ₂_guess::Real=0.02,gₛ₂_lim_lo::Real=0.001,gₛ₂_lim_hi::Real=0.5,
    Nₘ_f_guess::Real=0.012,Nₘ_f_lim_lo::Real=0.001,Nₘ_f_lim_hi::Real=0.05)
    
    x0 = [gₛ₁_guess, gₛ₂_guess, Nₘ_f_guess]    
  
    lower = [gₛ₁_lim_lo, gₛ₂_lim_lo, Nₘ_f_lim_lo]   
    upper = [gₛ₁_lim_hi, gₛ₂_lim_hi, Nₘ_f_lim_hi] 
    
    res = BlackBoxOptim.bboptimize(x->Objective_fun(x,daylength,photo_kinetic,envfun,model); SearchRange =[(gₛ₁_lim_lo,gₛ₁_lim_hi),(gₛ₂_lim_lo,gₛ₂_lim_hi),(Nₘ_f_lim_lo,Nₘ_f_lim_hi)])
    gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt = BlackBoxOptim.best_candidate(res)

    #=
    df = Optim.OnceDifferentiable(x->Objective_fun(x,daylength,photo_kinetic,envfun,model),x0;autodiff = :forward)   
    
    inner_optimizer = Optim.BFGS(linesearch = Optim.LineSearches.BackTracking())
    opt_trait = Optim.optimize(df, lower, upper, x0, Optim.Fminbox(inner_optimizer))
  
    Optim.converged(opt_trait)||error("No optimal traits could be found") 
      
    
    gₛ₁_opt = opt_trait.minimizer[1]
    gₛ₂_opt = opt_trait.minimizer[2]
    Nₘ_f_opt = opt_trait.minimizer[3]    
    =#
    
    return gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt
end

#initiate photo kinetic and environmental variables at time t (time after sunrise. Time in seconds).
function Init_weather_par!(t::Real,model::CCPHStruct,photo_kinetic::PhotoKineticRates,envfun::EnvironmentFunStruct)
    EnvironmentStruct!(t,model.env,envfun)
    PhotoPar!(model.photopar,photo_kinetic,model.env.Tₐ)
    return nothing
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
        
    return growthlength,step_length
end