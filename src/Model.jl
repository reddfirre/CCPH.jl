#Function for calculating Jmax at optimal (Jmaxₒₚₜ) temperature from Nitrogen per leaf area (Nₐ)
function Calc_Jmax(Nₐ::T,a::T,b::T) where {T<:Float64}     
    return max(a*Nₐ+b,0.0)
end
#Function for calcualting the seasonal peak Jmax from Nitrogen per leaf area (Nₐ)
function Calc_Jmax(Nₐ::T,a::T,b::T,b_Jmax::T) where {T<:Float64}     
    return b_Jmax*Calc_Jmax(Nₐ,a,b)
end
#Function for calcualting Jmax from Nitrogen per leaf area (Nₐ)
function Calc_Jmax(Nₐ::T,a::T,b::T,b_Jmax::T,xₜ::T) where {T<:Float64}     
    return xₜ*Calc_Jmax(Nₐ,a,b,b_Jmax) 
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
Calc_gₜ(gₛ::Float64,model::CCPHStruct) = gₛ*model.treepar.r_gₛ
Calc_gₜ(gₛ::Float64,treepar::TreePar) = gₛ*treepar.r_gₛ

#Calculate per sapwood mass nitrogen concentration
Calc_Nₘ_w(Nₘ_f::Float64,model::CCPHStruct) = model.treepar.rW*Nₘ_f
#Calcualte per fine roots mass nitrogen concentration
Calc_Nₘ_r(Nₘ_f::Float64,model::CCPHStruct) = model.treepar.rR*Nₘ_f      
#Calcualte per leaf area nitrogen concentration
Calc_Nₐ(Nₘ_f::Float64,model::CCPHStruct) = model.treepar.LMA*Nₘ_f

#Calculate total maintenance respiration
function Calc_Rₘ(Nₘ_f::T,Nₘ_w::T,Nₘ_r::T,Wf::T,Ww::T,Wr::T,model::CCPHStruct) where {T<:Float64}
    Rₘ = model.treepar.rₘ*(Nₘ_f*Wf+Nₘ_w*Ww+Nₘ_r*Wr)
    return Rₘ
end

Calc_NPP(P::T,Rₘ::T,model::CCPHStruct) where {T<:Float64} = model.treepar.y*(P-Rₘ)

#Foliage and fine root time depednent senescence
Calc_S_fr(Wf::T,Wr::T,model::CCPHStruct) where {T<:Float64} = Wf/model.treepar.Tf+Wr/model.treepar.Tr

#Tree height time derivative
function Calc_dH(NPP::T,S::T,Wf::T,Ww::T,Wr::T,H::T,Hs::T,model::CCPHStruct) where {T<:Float64}
    dH = (NPP-S)/( model.treepar.β₁*Ww/(model.treepar.β₁*H+model.treepar.β₂*Hs)+
    model.treepar.z*Wf/(H-Hs)+
    model.treepar.z*Ww/(H-Hs)+
    model.treepar.z*Wr/(H-Hs))

    return dH
end

#Function for calcualting fine root to sapwood area ratio based on carbon and nitrogen constraints
function Find_αr(Nₘ_f::T,Nₘ_w::T,Nₘ_r::T,P::T,model::CCPHStruct) where {T<:Float64}
    γC0 = model.treepar.y*(P-model.treepar.rₘ*(Nₘ_f*model.treesize.Wf+Nₘ_w*model.treesize.Ww))-model.treesize.Wf/model.treepar.Tf
    γC1 = model.treepar.y*model.treepar.rₘ*Nₘ_r*model.treesize.As+model.treesize.As/model.treepar.Tr
    γC2 = model.treepar.β₁*model.treesize.Ww/(model.treepar.β₁*model.treesize.H+model.treepar.β₂*model.treesize.Hs)+
        model.treepar.z*model.treesize.Wf/(model.treesize.H-model.treesize.Hs)+model.treepar.z*model.treesize.Ww/(model.treesize.H-model.treesize.Hs)
    γC3 = model.treepar.z*model.treesize.As/(model.treesize.H-model.treesize.Hs)   

    γN0 = Nₘ_f*model.treesize.Wf/model.treepar.Tf
    γN1 = Nₘ_r*model.treesize.As/model.treepar.Tr
    γN2 = Nₘ_w*model.treepar.β₁*model.treesize.Ww/(model.treepar.β₁*model.treesize.H+model.treepar.β₂*model.treesize.Hs)+
    Nₘ_f*model.treepar.z*model.treesize.Wf/(model.treesize.H-model.treesize.Hs)+Nₘ_w*model.treepar.z*model.treesize.Ww/(model.treesize.H-model.treesize.Hs)
    γN3 = Nₘ_r*model.treepar.z*model.treesize.As/(model.treesize.H-model.treesize.Hs) 

    γU0 = model.treepar.Nₛ*model.treesize.As*model.treepar.Kr
    γU1 = model.treepar.Kr
    γU2 = model.treesize.As*model.treesize.N

    θ₀ = (γC0*γN2+γN0*γC2)*γU1
    θ₁ = (γC0*γN3+γN1*γC2+γN0*γC3-γC1*γN2)*γU1+(γC0*γN2+γN0*γC2)*γU2-γC2*γU0
    θ₂ = (γC0*γN3+γN1*γC2+γN0*γC3-γC1*γN2)*γU2+(γN1*γC3-γC1*γN3)*γU1-γC3*γU0
    θ₃ = (γN1*γC3-γC1*γN3)*γU2

    if isnan(θ₀)||isnan(θ₁)||isnan(θ₂)||isnan(θ₃)
        error("Polynomial coefficients contains NaN")
    end
    if isinf(θ₀)||isinf(θ₁)||isinf(θ₂)||isinf(θ₃)
        error("Polynomial coefficients contains ∞")
    end
    
    αr_vec = solvecubic(θ₃, θ₂, θ₁, θ₀)     
    filter!(x->isnan(x)==false&&isinf(x)==false&&isreal(x)&&x>0.0,αr_vec)   

    isempty(αr_vec)==false||error("No feasible αr")    

    return αr_vec
end

#Calcualte gain (nett carbon gain minus fine root growth, performance measure)
function Gain_fun(αr::T,K_cost::T,Nₘ_f::T,Nₘ_w::T,Nₘ_r::T,P::T,model::CCPHStruct) where {T<:Float64}
    #Fine root mass
    Wr = model.treesize.As*αr
    
    #Calculate total maintenance respiration    
    R_m = Calc_Rₘ(Nₘ_f,Nₘ_w,Nₘ_r,model.treesize.Wf,model.treesize.Ww,Wr,model)

    #Net primary production   
    NPP = Calc_NPP(P,R_m,model)

    #Foliage and fine root time depednent senescence    
    S = Calc_S_fr(model.treesize.Wf,Wr,model)

    #Tree height time derivative
    dH = Calc_dH(NPP,S,model.treesize.Wf,model.treesize.Ww,Wr,model.treesize.H,model.treesize.Hs,model)

    #Fine root growth
    Gr = model.treepar.z*Wr/(model.treesize.H-model.treesize.Hs)*dH+Wr/model.treepar.Tr

    return (NPP-Gr)*K_cost^model.hydPar.i
end

#Run the Coupled Canopy Photosynthesis and Hydraulics model
function CCPH_run(gₛ::T,Nₘ_f::T,growthlength::T,model::CCPHStruct) where {T<:Float64}    
    #Calculate per sapwood mass nitrogen concentration    
    Nₘ_w = Calc_Nₘ_w(Nₘ_f,model)
    #Calcualte per fine roots mass nitrogen concentration       
    Nₘ_r = Calc_Nₘ_r(Nₘ_f,model)
    #Calcualte per leaf area nitrogen concentration    
    Nₐ = Calc_Nₐ(Nₘ_f,model)
    #Calculate Jmax
    Jmax = Calc_Jmax(Nₐ,model.treepar.a_Jmax,model.treepar.b_Jmax,model.photopar.b_Jmax,model.treepar.Xₜ)    
    #Quantum yield
    model.photopar.α = Calc_α(model.treepar.Xₜ,model.treepar.α_max)
    #Irradiance incident on a leaf at canopy top
    Iᵢ = Calc_Iᵢ(model.env.I₀,model)
    #Calculate LAI 
    LAI = Calc_LAI(model)
    #calcualte total conductance
    gₜ = Calc_gₜ(gₛ,model)
    #Calculate per tree carbon assimilation
    P = GPP(gₜ,Iᵢ,Jmax,LAI,growthlength,model)   
    #---Calculate cost factor of hydraulic failure---  
    K_cost, Kₓₗ, ψ_c = Calc_K_cost(gₛ,model)     
    #---Calculate posible fine root to sapwood area ratio, αr---    
    αr_vec = Find_αr(Nₘ_f,Nₘ_w,Nₘ_r,P,model)
    #Calcualte performance measure
    Gain_vec = [Gain_fun(αr,K_cost,Nₘ_f,Nₘ_w,Nₘ_r,P,model) for αr in αr_vec]

    #Find feasible αr which maximizes performance 
    filter!(x->isnan(x)==false&&isnan(x)==false&&isreal(x),Gain_vec)
    isempty(Gain_vec)==false||error("No feasible Gain")
    αr_val = αr_vec[argmax(Gain_vec)]

    modeloutput = CCPHOutput(P,αr_val,ψ_c,Kₓₗ,K_cost,maximum(Gain_vec))

    return modeloutput
end

function Trait_objective_fun(x::Array{T,1},growthlength::T,model::CCPHStruct) where {T<:Float64} 
    gₛ,Nₘ_f = x

    modeloutput = CCPH_run(gₛ,Nₘ_f,growthlength,model)

    return -modeloutput.Gain
end 

#Find optimal triats
function CCPHTraitmodel(growthlength::T,model::CCPHStruct;
    gₛ_guess::T=0.02,gₛ_lim_lo::T=0.001,gₛ_lim_hi::T=0.5,
    Nₘ_f_guess::T=0.012,Nₘ_f_lim_lo::T=0.001,Nₘ_f_lim_hi::T=0.1) where {T<:Float64}      
  
    x0 = [gₛ_guess, Nₘ_f_guess]    
  
    lower = [gₛ_lim_lo, Nₘ_f_lim_lo]   
    upper = [gₛ_lim_hi, Nₘ_f_lim_hi]
  
    inner_optimizer = Optim.BFGS(linesearch=Optim.LineSearches.BackTracking())   
    
    opt_trait = Optim.optimize(x->Trait_objective_fun(x,growthlength,model),
    lower,upper,x0,Optim.Fminbox(inner_optimizer))
  
    Optim.converged(opt_trait)||error("No optimal traits could be found")
    
    gₛ_opt = opt_trait.minimizer[1]
    Nₘ_f_opt = opt_trait.minimizer[2]
  
    return gₛ_opt,Nₘ_f_opt
end

#---This needs to be updated!---
#Differential equations describing the growth of the forest
function TreeStandDyn!(dy::Array{T,1},y::Array{T,1},model::CCPHStruct,
    t::T,growthlength::T,gₛ::T,Nₘ_f::T,αr::T) where {T<:Float64}  
    Wf,Ww,H,Hs,As,B,N,Wr = y

    #Calcualte per leaf area nitrogen concentration    
    Nₐ = Calc_Nₐ(Nₘ_f,model)
    #Calculate Jmax
    Jmax = Calc_Jmax(Nₐ,model.treepar.a_Jmax,model.treepar.b_Jmax,model.photopar.b_Jmax)
    #Irradiance incident on a leaf at canopy top
    Iᵢ = Calc_Iᵢ(model.env.I₀,model)
    #Calculate LAI 
    LAI = Calc_LAI(model)
    #calcualte total conductance
    gₜ = Calc_gₜ(gₛ,model)
    #Calculate GPP
    P = GPP(gₜ,Iᵢ,Jmax,LAI,growthlength,model)   

    #Calculate per sapwood mass nitrogen concentration    
    Nₘ_w = Calc_Nₘ_w(Nₘ_f,model)
    #Calcualte per fine roots mass nitrogen concentration       
    Nₘ_r = Calc_Nₘ_r(Nₘ_f,model)

    #Calculate total maintenance respiration    
    R_m = Calc_Rₘ(Nₘ_f,Nₘ_w,Nₘ_r,Wf,Ww,Wr,model)
    
    NPP = Calc_NPP(P,R_m,model)

    #Foliage and fine root time depednent senescence    
    S = Calc_S_fr(Wf,Wr,model)    

    #Calculate height growth  
    dH = Calc_dH(NPP,S,Wf,Ww,Wr,H,Hs,model)

    dN = 0 #Population density is constant 
    
    #---Calcualte leaf performance at crown base---
    Δ_leaf = Calc_Δ_leaf(gₜ,Iᵢ,LAI,growthlength,Nₘ_f,Jmax,model)

    dHs = dH*(Δ_leaf<0) #No crown rise when leaf performance >= 0

    dAs = model.treepar.z*As/(H-Hs)*(dH-dHs)
    dB = model.treepar.z*As/(H-Hs)*dH

    dWf = model.treepar.αf*dAs
    dWr = αr*dAs
    dWw = model.treepar.ρw*dAs*(model.treepar.β₁*H+model.treepar.β₂*Hs)+
    model.treepar.ρw*As*(model.treepar.β₁*dH+model.treepar.β₂*dHs)    
    
    dy[1] = dWf
    dy[2] = dWw
    dy[3] = dH
    dy[4] = dHs
    dy[5] = dAs  
    dy[6] = dB
    dy[7] = dN
    dy[8] = dWr
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

#---This needs to be updated!---
#Simulate growth with weather time series
function CCPHStandGrowth!(model::CCPHStruct,weatherts::WeatherTS,
    photo_kinetic::PhotoKineticRates,nstep::Integer;
    K_cost_crit::T=0.12,gₛ_guess::T=0.02,Nₘ_f_guess::T=0.012)::CCPHTS  where {T<:Float64}
    
    ccphts = CCPHTS() #Init Time series struct
    CCPHTS!(ccphts,model.treesize)      
       
    for i in 1:nstep
                                   
        #initiate weather parameters         
        growthlength,step_length = Init_weather_par!(i,model,weatherts,photo_kinetic)

        gₛ_crit = Calc_K_costᵢₙᵥ(K_cost_crit,model)         

        gₛ_opt,Nₘ_f_opt = CCPHTraitmodel(growthlength,model;
        gₛ_guess=gₛ_guess,Nₘ_f_guess=Nₘ_f_guess,gₛ_lim_hi=gₛ_crit)
          
        modeloutput = CCPH_run(gₛ_opt,Nₘ_f_opt,growthlength,model)    
      
        CCPHTS!(ccphts,modeloutput,gₛ_opt,Nₘ_f_opt,gₛ_crit)         
  
        y0 = [model.treesize.Wf,model.treesize.Ww,model.treesize.H,
        model.treesize.Hs,model.treesize.As,model.treesize.B,model.treesize.N,modeloutput.αr*model.treesize.As]        
        
        #Simulate one timestep
        prob = DifferentialEquations.ODEProblem((dy::Array{Float64,1},y::Array{Float64,1},P_in::CCPHStruct,t::Float64)->
        TreeStandDyn!(dy,y,P_in,t,growthlength,gₛ_opt,Nₘ_f_opt,modeloutput.αr),y0,(0.0,step_length),model)
        sol = DifferentialEquations.solve(prob)     
  
        #Update current size        
        model.treesize.Wf,model.treesize.Ww,model.treesize.H,
        model.treesize.Hs,model.treesize.As,model.treesize.B,model.treesize.N,Wr = sol(step_length)          
  
        #Update time series
        CCPHTS!(ccphts,model.treesize)         
    end    

    #initiate weather parameters
    growthlength,step_length = Init_weather_par!(nstep+1,model,weatherts,photo_kinetic)

    gₛ_crit = Calc_K_costᵢₙᵥ(K_cost_crit,model)
    
    gₛ_opt,Nₘ_f_opt = CCPHTraitmodel(growthlength,model;
    gₛ_guess=gₛ_guess,Nₘ_f_guess=Nₘ_f_guess,gₛ_lim_hi=gₛ_crit)    
    
    modeloutput = CCPH_run(gₛ_opt,Nₘ_f_opt,growthlength,model)        
    
    CCPHTS!(ccphts,modeloutput,gₛ_opt,Nₘ_f_opt,gₛ_crit)

    return ccphts
end