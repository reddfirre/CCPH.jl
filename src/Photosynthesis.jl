#Rate of electron transport
function Calc_J(Iᵢ::T,Jₘₐₓ::S,α::T,θ::T) where {S<:Real,T<:Float64} 
    J = (α*Iᵢ+Jₘₐₓ-sqrt(α^2*Iᵢ^2+2*α*Iᵢ*Jₘₐₓ*(1-2*θ)+Jₘₐₓ^2))/(2*θ)

    return J
end
#Function for calculating intercellular carbon dioxide concentration  
function calc_opt_cᵢ(gₜ::S, 
    cₐ::T,
    J::S,
    Γ::T,
    Pₜₒₜ::T) where {T<:Real,S<:Real}
    gₜₚ = gₜ/Pₜₒₜ
    a₀ = J/(4*gₜₚ)*Γ+2*cₐ*Γ
    a₁ = cₐ-2*Γ-J/(4*gₜₚ)
    cᵢ = a₁/2+sqrt(a₁^2/4+a₀)
    Γ<cᵢ<cₐ||error("calc_opt_cᵢ: Γ<cᵢ<cₐ failed")
    return cᵢ
end
function calc_opt_cᵢ(gₜ::T,    
    J::T,
    photo::PhotoPar,
    env::EnvironmentStruct) where {T<:Real}

    Pₜₒₜ,cₐ,Γ = env.P,env.Cₐ,photo.Γ    

    cᵢ = calc_opt_cᵢ(gₜ,cₐ,J,Γ,Pₜₒₜ)
    return cᵢ
end
#Calculate assimilation rate (mol C m⁻² leaf area s⁻¹)
function calc_Assimilation(gₜ::S,cᵢ::S,Pₜₒₜ::T,cₐ::T) where {T<:Real,S<:Real}
    A = gₜ*(cₐ-cᵢ)/Pₜₒₜ
    return A
end
#Calculate Vcmax
function calc_Vcmax(J::S,
    cᵢ::S,
    K::T,
    Γ::T) where {T<:Real,S<:Real}
    Vcmax = J*(cᵢ+K)/(4*(cᵢ+2*Γ))
    return Vcmax
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

#Calculate per tree canopy gross primary production (kg C year⁻¹ tree⁻¹)
GPP(A::S,LAI::T,growthlength::T,model::CCPHStruct) where {S<:Real,T<:Real} = 
model.cons.M_C*A*growthlength*(1-exp(-model.treepar.k*LAI))/(model.treesize.N*model.treepar.k)

#=
function GPP(gₛ::S,Nₘ_f::S,growthlength::T,model::CCPHStruct) where {S<:Real,T<:Float64}
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
    #Calculate C assimilation
    A = calc_Assimilation(gₜ,cᵢ,model.env.P,model.env.Cₐ)
    #Calculate per tree carbon assimilation
    P =GPP(A,LAI,growthlength,model)
    return P
end
=#