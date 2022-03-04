#vulnerability curve
Pfun(ψ::T,ψ₅₀::T,b::T) where {T<:Float64} = (1/2)^(ψ/ψ₅₀)^b

#Invers vulnerability curve
Pfunᵢₙᵥ(Pval::T,ψ₅₀::T,b::T) where {T<:Float64} = ψ₅₀*(log(Pval)/log(0.5))^(1/b)

#Effective saturation 
CalcSₑ(θₛ::T;θₛₐₜ::T=0.41,θᵣ::T=0.006) where {T<:Float64} = (θₛ-θᵣ)/(θₛₐₜ-θᵣ)

#Soil water content (volumetric) to soil water potential (MPa) (Water retention curve) 
function θₛ2ψₛ(θₛ::T;θₛₐₜ::T=0.41,θᵣ::T=0.006,λ::T=1.0,ψₐ::T=-0.097706) where {T<:Float64}
    Sₑ = CalcSₑ(θₛ;θₛₐₜ=θₛₐₜ,θᵣ=θᵣ)

    return ψₐ*Sₑ^(-1/λ)
end    

#Relative Soil conductance (ratio between actual and maximal)
function Re_Kₛᵣfun(θₛ::T;θₛₐₜ::T=0.41,θᵣ::T=0.006,p::T=4.66) where {T<:Float64}
    Sₑ = CalcSₑ(θₛ;θₛₐₜ=θₛₐₜ,θᵣ=θᵣ)
    return Sₑ^p
end

#Calculate canopy conductance
function Calc_K_cost(gₛ::T,model::CCPHStruct;limit_up::T=1.0,limit_lo::T=0.12) where {T<:Float64}
    ψ₅₀,b,Kₓₗ₀,g,ρ_H2O,θₛ = model.hydPar.ψ₅₀,model.hydPar.b,model.hydPar.Kₓₗ₀,model.cons.g,model.cons.ρ_H2O,model.env.θₛ

    #Calculate tree height
    H = model.treesize.H
    #Calculate transpiration    
    E = 1.6*gₛ*model.env.VPD/model.env.P
    #Caluclate soil water potential
    ψₛ = θₛ2ψₛ(θₛ)      

    #Soil water potential adjusted for gravitational pressure (MPa)
    ψₛ_g = ψₛ-H*ρ_H2O*g*10^-6  

    #Calculate canopy conductance
    try
        K_cost = bisection(x->x-Pfun(ψₛ_g-E/(2*Kₓₗ₀*x),ψ₅₀,b),limit_lo, limit_up) #lower than 0.12 resultts in hydraulic failure
        Kₓₗ = Kₓₗ₀*K_cost
        ψ_c = ψₛ_g-E/Kₓₗ       
        return K_cost, Kₓₗ, ψ_c 
    catch
        error("Could not find a feasable canopy conductance")
    end    
end

#Calculate stomatal (gₛ) for a given K_cost (Pval)
function Calc_K_costᵢₙᵥ(Pval::Float64,model::CCPHStruct)
    ψ₅₀,b,Kₓₗ₀,g,ρ_H2O,θₛ = model.hydPar.ψ₅₀,model.hydPar.b,model.hydPar.Kₓₗ₀,model.cons.g,model.cons.ρ_H2O,model.env.θₛ

    #Calc target ψ
    ψ_target = Pfunᵢₙᵥ(Pval,ψ₅₀,b)
    #Calc taget Kₓₗ
    Kₓₗ = Kₓₗ₀*Pval
    #Calculate tree height
    H = model.treesize.H    
    #Caluclate soil water potential
    ψₛ = θₛ2ψₛ(θₛ)    
    #Soil water potential adjusted for gravitational pressure (MPa)
    ψₛ_g = ψₛ-H*ρ_H2O*g*10^-6     
    #Calculate target E
    E_target = -(ψ_target-ψₛ_g)*2*Kₓₗ
    #Calcualte target gₛ
    gₛ_target = E_target*model.env.P/(1.6*model.env.VPD)

    return gₛ_target
end