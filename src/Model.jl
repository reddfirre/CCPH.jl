#Function for calcualting Jmax from Nitrogen per leaf area (Nₐ)
function Calc_Jmax(Nₐ::T,a::T,b::T) where {T<:Float64}     
    return max(a*Nₐ+b,0.0)
end

Calc_LAI(model::CCPHStruct) = model.treesize.Wf/model.treepar.LMA*model.treesize.N

#Calculate the irradiance incident on a leaf at irradiance I
Calc_Iᵢ(I::Float64,model::CCPHStruct) = I*model.treepar.k/(1-model.treepar.m)

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
    R_m = model.treepar.rₘ*(Nₘ_f*model.treesize.Wf+Nₘ_w*model.treesize.Ww+Nₘ_r*Wr)

    #Net primary production
    NPP = model.treepar.y*(P-R_m)

    #Foliage and fine root time depednent senescence
    S = model.treesize.Wf/model.treepar.Tf+Wr/model.treepar.Tr

    #Tree height time derivative
    dH = (NPP-S)/( model.treepar.β₁*model.treesize.Ww/(model.treepar.β₁*model.treesize.H+model.treepar.β₂*model.treesize.Hs)+
    model.treepar.z*model.treesize.Wf/(model.treesize.H-model.treesize.Hs)+
    model.treepar.z*model.treesize.Ww/(model.treesize.H-model.treesize.Hs)+
    model.treepar.z*Wr/(model.treesize.H-model.treesize.Hs))

    #Fine root growth
    Gr = model.treepar.z*Wr/(model.treesize.H-model.treesize.Hs)*dH+Wr/model.treepar.Tr

    return (NPP-Gr)*K_cost^model.hydPar.i
end

function Simple_CCPH(gₛ::T,Nₘ_f::T,growthlength::T,model::CCPHStruct) where {T<:Float64}    
    #Calculate per sapwood mass nitrogen concentration
    Nₘ_w = model.treepar.rW*Nₘ_f
    #Calcualte per fine roots mass nitrogen concentration
    Nₘ_r = model.treepar.rR*Nₘ_f      
    #Calcualte per leaf area nitrogen concentration
    Nₐ = model.treepar.LMA*Nₘ_f
    #Calculate Jmax
    Jmax = Calc_Jmax(Nₐ,model.treepar.a_Jmax,model.treepar.b_Jmax)
    #Irradiance incident on a leaf at canopy top
    Iᵢ = Calc_Iᵢ(model.env.I₀,model)
    #Calculate LAI 
    LAI = Calc_LAI(model)
    #calcualte totel conductance
    gₜ = gₛ*model.treepar.r_gₛ
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