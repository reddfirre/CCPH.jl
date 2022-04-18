#Rate of electron transport
function Calc_J(Iᵢ::T,Jₘₐₓ::S,α::T,θ::T) where {S<:Real,T<:Float64} 
    J = (α*Iᵢ+Jₘₐₓ-sqrt(α^2*Iᵢ^2+2*α*Iᵢ*Jₘₐₓ*(1-2*θ)+Jₘₐₓ^2))/(2*θ)

    return J
end

#Photosynthesis model for given total conductance (stomatal+mesophyll), gₜ, irradiance incident on the leaf (Iᵢ), 
#and Jₘₐₓ
function Farquhar(gₜ::S, Iᵢ::T, Jₘₐₓ::S,photo::PhotoPar,env::EnvironmentStruct) where {S<:Real,T<:Float64}   
    Pₜₒₜ,cₐ = env.P,env.Cₐ
    
    gₜₚ = gₜ/Pₜₒₜ
    
    J = Calc_J(Iᵢ,Jₘₐₓ,photo.α,photo.θ)
    Aⱼ = J/4

    #Intercellular carbon dioxide concentration
    cᵢ = (gₜₚ*(cₐ-2*photo.Γ)-Aⱼ*(1-photo.b_r)+sqrt((gₜₚ*(cₐ-2*photo.Γ)-Aⱼ*(1-photo.b_r))^2+4*Aⱼ*gₜₚ*(photo.Γ+photo.b_r*photo.K)+8*gₜₚ^2*cₐ*photo.Γ))/(2*gₜₚ) 
     
    #Carbon assimilation rate
    A = gₜₚ*(cₐ-cᵢ)  

    return A, cᵢ
end

#Calculate the C assimilation (kg C year⁻¹ m⁻² leaf area)
function C_assimilation(gₜ::S,Iᵢ::T,Jmax::S,growthlength::T,model::CCPHStruct) where {S<:Real,T<:Float64}
    A = model.cons.M_C*Farquhar(gₜ,Iᵢ,Jmax,model.photopar,model.env)[1]*growthlength

    return A
end
#Calculate per tree canopy gross primary production (kg C year⁻¹ tree⁻¹)
GPP(gₜ::S,Iᵢ::T,Jmax::S,LAI::T,growthlength::T,model::CCPHStruct) where {S<:Real,T<:Float64} = 
C_assimilation(gₜ,Iᵢ,Jmax,growthlength,model)*(1-exp(-model.treepar.k*LAI))/(model.treesize.N*model.treepar.k)
function GPP(gₛ::S,Nₘ_f::S,growthlength::T,model::CCPHStruct) where {S<:Real,T<:Float64}
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
     return P
end

#Calcualte upper limit for GPP
function GPP_max(gₛ::T,growthlength::T,model::CCPHStruct) where {T<:Float64}
    k,N = model.treepar.k,model.treesize.N
    #Calculate LAI 
    LAI = Calc_LAI(model)
    #calcualte total conductance
    gₜ = Calc_gₜ(gₛ,model)    
    cᵢ_min = (model.photopar.Γ+model.photopar.b_r*model.photopar.K)/(1-model.photopar.b_r)
    GPP = gₜ*(model.env.Cₐ-cᵢ_min)*(1-exp(-k*LAI))/(N*k)*model.cons.M_C*growthlength/model.env.P
    return GPP
end

#---This needs to be updated!---
#Calcualte leaf performance at crown base
function Calc_Δ_leaf(gₜ::T,Iᵢ::T,LAI::T,growthlength::T,Nₘ_f::T,Jmax::T,model::CCPHStruct) where {T<:Float64} 
    Iᵢ_b = Iᵢ*exp(-model.treepar.k*LAI)
    Jmax_b = Jmax*exp(-model.treepar.k*LAI)
    A_b = C_assimilation(gₜ,Iᵢ_b,Jmax_b,growthlength,model)  
    Δ_leaf = model.treepar.y*(A_b-model.treepar.rₘ*Nₘ_f*model.treepar.LMA)-model.treepar.LMA/model.treepar.Tf #Bottom leaf performance
    return Δ_leaf
end