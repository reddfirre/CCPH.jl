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
#Calculate instantaneous leaf assimilation rate (mol C m⁻² leaf area s⁻¹)
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

#Calculate per tree canopy gross primary production (kg C year⁻¹ tree⁻¹)
GPP(A::S,LAI::T,growthlength::T,model::CCPHStruct) where {S<:Real,T<:Real} = 
model.cons.M_C*A*growthlength*(1-exp(-model.treepar.k*LAI))/(model.treesize.N*model.treepar.k)