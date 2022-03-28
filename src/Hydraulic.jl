#vulnerability curve
Pfun(ψ::T,ψ₅₀::T,b::T) where {T<:Float64} = (1/2)^(ψ/ψ₅₀)^b

#Invers vulnerability curve
Pfunᵢₙᵥ(Pval::T,ψ₅₀::T,b::T) where {T<:Float64} = ψ₅₀*(log(Pval)/log(0.5))^(1/b)

#The integral of Pfun from ψ to -∞
Pintlim(ψ::T,ψ₅₀::T,b::T) where {T<:Float64} = CCPH.Pfun(ψ,ψ₅₀,b)*ψ+
abs(ψ₅₀)/log(2)^(1/b)*SpecialFunctions.gamma((1+b)/b,log(2)*(ψ/ψ₅₀)^b)

#The integral of Pfun from the pre-dawn canopy water potential ψ_cm to the canopy water potential ψ_c 
Pint(ψ_c::T,ψ_cm::T,ψ₅₀::T,b::T) where {T<:Float64} = Pintlim(ψ_cm,ψ₅₀,b)-Pintlim(ψ_c,ψ₅₀,b)

#Fix-point equation for calculating the xylem conductance Kₓₗ=Kₓₗ₀*x
Ptarget(x::T,Kₓₗ₀::T,E::T,ψ_cm::T,ψ₅₀::T,b::T) where {T<:Float64} = Pint(ψ_cm-E/(Kₓₗ₀*x),ψ_cm,ψ₅₀,b)/(E/(Kₓₗ₀*x))

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
function Calc_K_cost(gₛ::T,model::CCPHStruct) where {T<:Float64}
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
        K_cost = fixpoint(x->Ptarget(x,Kₓₗ₀,E,ψₛ_g,ψ₅₀,b),1.0)
        Kₓₗ = Kₓₗ₀*K_cost
        ψ_c = ψₛ_g-E/Kₓₗ 
        return K_cost, Kₓₗ, ψ_c     
    catch
        error("K_cost failed: gₛ=$(gₛ), E=$(E), Kₓₗ₀=$(Kₓₗ₀), ψₛ_g =$(ψₛ_g)")
    end      
end

#Calculate the stomatal conductance (gₛ/E) such that Pfun(ψ_target,ψ₅₀,b)=Pval and E = Pint(ψ_target,ψₛ_g,ψ₅₀,b)
function Calc_K_costᵢₙᵥ(Pval::Float64,model::CCPHStruct)
    ψ₅₀,b,Kₓₗ₀,g,ρ_H2O,θₛ = model.hydPar.ψ₅₀,model.hydPar.b,model.hydPar.Kₓₗ₀,model.cons.g,model.cons.ρ_H2O,model.env.θₛ

    #Calc target ψ
    ψ_target = Pfunᵢₙᵥ(Pval,ψ₅₀,b)    
    #Calculate tree height
    H = model.treesize.H    
    #Caluclate soil water potential
    ψₛ = θₛ2ψₛ(θₛ)    
    #Soil water potential adjusted for gravitational pressure (MPa)
    ψₛ_g = ψₛ-H*ρ_H2O*g*10^-6     
    #Calculate target E
    E_target = Kₓₗ₀*Pint(ψ_target,ψₛ_g,ψ₅₀,b)
    #Calcualte target gₛ
    gₛ_target = E_target*model.env.P/(1.6*model.env.VPD)

    return gₛ_target
end