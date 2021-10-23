#Photosynthesis model for given total conductance (stomatal+mesophyll), gₜ, irradiance incident on the leaf (Iᵢ), 
#and Jₘₐₓ
function Farquhar(gₜ::T, Iᵢ::T, Jₘₐₓ::T,photo::PhotoPar,env::EnvironmentStruct) where {T<:Float64}   
    Pₜₒₜ,cₐ = env.P,env.Cₐ
    
    gₜₚ = gₜ/Pₜₒₜ
    
    J = (photo.α*Iᵢ+Jₘₐₓ-sqrt(photo.α^2*Iᵢ^2+2*photo.α*Iᵢ*Jₘₐₓ*(1-2*photo.θ)+Jₘₐₓ^2))/(2*photo.θ)
    Aⱼ = J/4

    #Intercellular carbon dioxide concentration
    cᵢ = (gₜₚ*(cₐ-2*photo.Γ)-Aⱼ*(1-photo.b_r)+sqrt((gₜₚ*(cₐ-2*photo.Γ)-Aⱼ*(1-photo.b_r))^2+4*Aⱼ*gₜₚ*(photo.Γ+photo.b_r*photo.K)+8*gₜₚ^2*cₐ*photo.Γ))/(2*gₜₚ) 
     
    #Carbon assimilation rate
    A = gₜₚ*(cₐ-cᵢ)  

    return A, cᵢ
end

#Calculate per tree canopy gross primary production
GPP(gₜ::T,Iᵢ::T,Jmax::T,LAI::T,growthlength::T,model::CCPHStruct) where {T<:Float64} = model.cons.M_C*Farquhar(gₜ,Iᵢ,Jmax,model.photopar,model.env)[1]*
    growthlength*(1-exp(-model.treepar.k*LAI))/(model.treesize.N*model.treepar.k)