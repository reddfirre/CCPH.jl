#Function for calculating Jmax at optimal (Jmaxₒₚₜ) temperature from Nitrogen per leaf area (Nₐ)
function Calc_Jmax(Nₐ::S,a::T,b::T) where {S<:Real,T<:Float64}     
    return max(a*Nₐ+b,0.0)
end
#Function for calcualting the seasonal peak Jmax from Nitrogen per leaf area (Nₐ)
function Calc_Jmax(Nₐ::S,a::T,b::T,b_Jmax::T) where {S<:Real,T<:Float64}     
    return b_Jmax*Calc_Jmax(Nₐ,a,b)
end
#Function for calcualting Jmax from Nitrogen per leaf area (Nₐ)
function Calc_Jmax(Nₐ::S,a::T,b::T,b_Jmax::T,xₜ::T) where {S<:Real,T<:Float64}     
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
Calc_gₜ(gₛ::Real,model::CCPHStruct) = gₛ*model.treepar.r_gₛ
Calc_gₜ(gₛ::Real,treepar::TreePar) = gₛ*treepar.r_gₛ

#Calculate per sapwood mass nitrogen concentration
Calc_Nₘ_w(Nₘ_f::Real,model::CCPHStruct) = model.treepar.rW*Nₘ_f
#Calcualte per fine roots mass nitrogen concentration
Calc_Nₘ_r(Nₘ_f::Real,model::CCPHStruct) = model.treepar.rR*Nₘ_f      
#Calcualte per leaf area nitrogen concentration
Calc_Nₐ(Nₘ_f::Real,model::CCPHStruct) = model.treepar.LMA*Nₘ_f

#Calcualte the ration between Nₘ and Nₘ_f
function Calc_Δₘ(H::T,Hₛ::T,Wf::T,Ww::T,Wr::T,r_w::T,r_r::T,β₁::T,β₂::T,z::T) where {T<:Float64}
    Δ = z/(H-Hₛ)
    Δw = β₁*Ww/(β₁*H+β₂*Hₛ)

    Δₘ = (r_w*Δw+Δ*(Wf+r_w*Ww+r_r*Wr))/(Δw+Δ*(Wf+Ww+Wr))
    return Δₘ
end
function Calc_Δₘ(Wr::T,model::CCPHStruct) where {T<:Float64}
    r_w = model.treepar.rW
    r_r = model.treepar.rR
    H,Hₛ,Wf,Ww = model.treesize.H,model.treesize.Hs,model.treesize.Wf,model.treesize.Ww
    β₁,β₂,z = model.treepar.β₁,model.treepar.β₂,model.treepar.z

    Δₘ = Calc_Δₘ(H,Hₛ,Wf,Ww,Wr,r_w,r_r,β₁,β₂,z)
    return Δₘ
end
function Calc_Δₘ(model::CCPHStruct)   

    Δₘ = Calc_Δₘ(0.0,model)
    return Δₘ
end

#Parameter values for drawing the growth constraint function
function Calc_Par(P::T,Nₘ_f::T,model::CCPHStruct) where {T<:Real} 
    H,Hₛ,Wf,Ww,N = model.treesize.H,model.treesize.Hs,model.treesize.Wf,model.treesize.Ww,model.treesize.N
    β₁,β₂,z,Nₛ = model.treepar.β₁,model.treepar.β₂,model.treepar.z,model.treepar.Nₛ
    Tᵣ,rₘ,y = model.treepar.Tr,model.treepar.rₘ,model.treepar.y
    Nₘ_w = Calc_Nₘ_w(Nₘ_f,model)
    NPPΔ = y*(P-rₘ*Nₘ_f*Wf-rₘ*Nₘ_w*Ww)
    H₁ = z/(H-Hₛ)
    H₂ = β₁*Ww/(β₁*H+β₂*Hₛ)
    return (NPPΔ,Nₘ_f,rₘ,Wf,Ww,H₁,H₂,Nₛ,N)
end

#Calculate total maintenance respiration
function Calc_Rₘ(Nₘ_f::S,Nₘ_w::S,Nₘ_r::S,Wf::T,Ww::T,Wr::S,model::CCPHStruct) where {S<:Real,T<:Float64}
    Rₘ = model.treepar.rₘ*(Nₘ_f*Wf+Nₘ_w*Ww+Nₘ_r*Wr)
    return Rₘ
end

Calc_NPP(P::T,Rₘ::T,model::CCPHStruct) where {T<:Real} = model.treepar.y*(P-Rₘ)

#Foliage and fine root time depednent senescence
Calc_S_fr(Wf::T,Wr::S,model::CCPHStruct) where {S<:Real,T<:Float64} = Wf/model.treepar.Tf+Wr/model.treepar.Tr

#Tree height time derivative
function Calc_dH(NPP::W,S::W,Wf::T,Ww::T,Wr::W,H::T,Hs::T,model::CCPHStruct) where {W<:Real,T<:Float64}
    dH = (NPP-S)/( model.treepar.β₁*Ww/(model.treepar.β₁*H+model.treepar.β₂*Hs)+
    model.treepar.z*Wf/(H-Hs)+
    model.treepar.z*Ww/(H-Hs)+
    model.treepar.z*Wr/(H-Hs))

    return dH
end

#Function for calcualting fine root to sapwood area ratio based on carbon and nitrogen constraints
function Find_αr(Nₘ_f::T,Nₘ_w::T,Nₘ_r::T,P::T,model::CCPHStruct) where {T<:Real}
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
        error("Find_αr: Polynomial coefficients contains NaN")
    end
    if isinf(θ₀)||isinf(θ₁)||isinf(θ₂)||isinf(θ₃)
        error("Find_αr: Polynomial coefficients contains ∞")
    end
    
    αr_vec = solvecubic(θ₃, θ₂, θ₁, θ₀)     
    filter!(x->isnan(x)==false&&isinf(x)==false&&isreal(x)&&x>0.0,αr_vec)   

    isempty(αr_vec)==false||error("Find_αr: No feasible αr: (NPPΔ,Nₘ_f,rₘ,Wf,Ww,H₁,H₂,Nₛ,N)=$(Calc_Par(P,Nₘ_f,model))")    

    return αr_vec
end

#Calcualte gain (nett carbon gain minus fine root growth, performance measure)
function Gain_fun(αr::T,K_cost::T,Nₘ_f::T,Nₘ_w::T,Nₘ_r::T,P::T,model::CCPHStruct) where {T<:Real}
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
function CCPH_run(gₛ::S,Nₘ_f::S,growthlength::T,model::CCPHStruct) where {S<:Real,T<:Float64}    
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

#Nₘ_f<Nₘ_f_pos guarantees solution to Find_αr 
function Calc_Nₘ_f_pos(model::CCPHStruct)
    Δₘ₀ = Calc_Δₘ(model)
    Wf,Ww = model.treesize.Wf,model.treesize.Ww
    Tf,rₘ,y,rW = model.treepar.Tf,model.treepar.rₘ,model.treepar.y,model.treepar.rW

    Nₘ_f_pos = Wf*(1-Δₘ₀)/(Tf*y*rₘ*Δₘ₀*(Wf+rW*Ww))
    return Nₘ_f_pos
end

#Calcualte Δ. Δ>0 guarantees a solution for Find_αr
function Calc_Δ(x::S,growthlength::T,model::CCPHStruct) where {S<:AbstractArray,T<:Float64}
    gₛ,Nₘ_f = x
    Wf,Ww = model.treesize.Wf,model.treesize.Ww
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
    Δₘ₀ = Calc_Δₘ(model) 
    Δ = model.treepar.y*(P-model.treepar.rₘ*(Nₘ_f*Wf+Nₘ_w*Ww))*Δₘ₀+Wf*(1-Δₘ₀)/model.treepar.Tf
    
    return Δ
end
#Calculate trait optimization constraint
con!(c::S,x::S,growthlength::T,model::CCPHStruct) where {S<:AbstractArray,T<:Float64} = 
(c[1] = Calc_Δ(x,growthlength,model);c)
#Calculate trait optimization constraint Jacobian
function con_J!(J::AbstractMatrix,x::S,growthlength::T,model::CCPHStruct) where {S<:AbstractArray,T<:Float64}
    J[1,1],J[1,2] = ForwardDiff.gradient(y->Calc_Δ(y,growthlength,model),x)
    J
end
#Calculate trait optimization constraint Hessian
function con_H!(H::AbstractMatrix,x::S,λ::S,growthlength::T,model::CCPHStruct) where {S<:AbstractArray,T<:Float64}
    H_temp = ForwardDiff.hessian(y->Calc_Δ(y,growthlength,model),x)
    H[1,1] += λ[1]*H_temp[1,1]
    H[1,2] += λ[1]*H_temp[1,2]
    H[2,1] += λ[1]*H_temp[2,1]
    H[2,2] += λ[1]*H_temp[2,2]
    H    
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
  
    Nₘ_f_pos = Calc_Nₘ_f_pos(model)

    x0 = [gₛ_guess, min(Nₘ_f_guess,Nₘ_f_pos)]    
  
    lower = [gₛ_lim_lo, min(Nₘ_f_lim_lo,Nₘ_f_pos*0.8)]   
    upper = [min(gₛ_lim_hi,0.5), Nₘ_f_lim_hi] 
    lc = [0.0]; uc = [Inf] 

    df = Optim.TwiceDifferentiable(x->Trait_objective_fun(x,growthlength,model),x0;autodiff = :forward)   
    
    dfc = Optim.TwiceDifferentiableConstraints((c,x)->con!(c,x,growthlength,model),
    (J,x)->con_J!(J,x,growthlength,model),
    (H,x,λ)->con_H!(H,x,λ,growthlength,model),
    lower,upper,lc,uc)  

    option = Optim.Options(allow_f_increases = true, successive_f_tol = 2,f_tol=10^-6,time_limit=30.0)

    opt_trait = Optim.optimize(df, dfc, x0, Optim.IPNewton(),option)
  
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