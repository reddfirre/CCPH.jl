#vulnerability curve
Pfun(ψ::S,ψ₅₀::T,b::T) where {S<:Real,T<:Float64} = (1/2)^(ψ/ψ₅₀)^b

#Invers vulnerability curve
Pfunᵢₙᵥ(Pval::T,ψ₅₀::T,b::T) where {T<:Float64} = ψ₅₀*(log(Pval)/log(0.5))^(1/b)

#The integral of Pfun from ψ to -∞
Pintlim(ψ::S,ψ₅₀::T,b::T) where {S<:Real,T<:Float64} = CCPH.Pfun(ψ,ψ₅₀,b)*ψ+
abs(ψ₅₀)/log(2)^(1/b)*SpecialFunctions.gamma((1+b)/b,log(2)*(ψ/ψ₅₀)^b)

function Pintlim(ψ::S,ψ₅₀::T,b::T,Γₐ::T) where {S<:Real,T<:Real} 
    p,q = SpecialFunctions.gamma_inc((1+b)/b,log(2)*(ψ/ψ₅₀)^b,2)
    return CCPH.Pfun(ψ,ψ₅₀,b)*ψ+abs(ψ₅₀)/log(2)^(1/b)*Γₐ*q
end

#The integral of Pfun from the pre-dawn canopy water potential ψ_cm to the canopy water potential ψ_c 
#Pint(ψ_c::S,ψ_cm::T,ψ₅₀::T,b::T) where {S<:Real,T<:Float64} = Pintlim(ψ_cm,ψ₅₀,b)-Pintlim(ψ_c,ψ₅₀,b)

function Pint(ψ_c::S,ψ_cm::T,ψ₅₀::T,b::T)  where {S<:Real,T<:Real}
    ψ = (ψ_c+ψ_cm)/2
    Δψ = (ψ_cm-ψ_c)
    return Δψ*(CCPH.Pfun(ψ_c,ψ₅₀,b)+4*CCPH.Pfun(ψ,ψ₅₀,b)+CCPH.Pfun(ψ_cm,ψ₅₀,b))/6
end

Pint(ψ_c::S,ψ_cm::T,ψ₅₀::T,b::T,Γₐ::T) where {S<:Real,T<:Real} = Pintlim(ψ_cm,ψ₅₀,b,Γₐ)-Pintlim(ψ_c,ψ₅₀,b,Γₐ)

#Fix-point equation for calculating the xylem conductance Kₓₗ=Kₓₗ₀*x
Ptarget(x::W,Kₓₗ₀::T,E::S,ψ_cm::T,ψ₅₀::T,b::T) where {S<:Real,W<:Real,T<:Float64} = Pint(ψ_cm-E/(Kₓₗ₀*x),ψ_cm,ψ₅₀,b)/(E/(Kₓₗ₀*x))

Ptarget(x::W,Kₓₗ₀::T,E::S,ψ_cm::T,ψ₅₀::T,b::T,Γₐ::T) where {S<:Real,W<:Real,T<:Float64} = Pint(ψ_cm-E/(Kₓₗ₀*x),ψ_cm,ψ₅₀,b,Γₐ::T)/(E/(Kₓₗ₀*x))

#Find root to equation to find the xylem conductance Kₓₗ=Kₓₗ₀*x
P_zero(x::W,Kₓₗ₀::T,E::S,ψ_cm::T,ψ₅₀::T,b::T) where {S<:Real,W<:Real,T<:Real} = Pint(ψ_cm-E/(Kₓₗ₀*x),ψ_cm,ψ₅₀,b)*Kₓₗ₀/E-1

#Effective saturation 
calcSₑ(θₛ::T;θₛₐₜ::T=0.41,θᵣ::T=0.006) where {T<:Float64} = (θₛ-θᵣ)/(θₛₐₜ-θᵣ)

#Soil water content (volumetric) to soil water potential (MPa) (Water retention curve) 
function θₛ2ψₛ(θₛ::T;θₛₐₜ::T=0.41,θᵣ::T=0.006,λ::T=1.0,ψₐ::T=-0.097706) where {T<:Float64}
    Sₑ = calcSₑ(θₛ;θₛₐₜ=θₛₐₜ,θᵣ=θᵣ)

    return ψₐ*Sₑ^(-1/λ)
end    

#Calculate leaf transpiration E (mol H₂O m⁻² leaf area s⁻¹)
function calc_E(gₛ::S,
    VPD::T,
    P::T;
    cons::Constants=Constants()) where {T<:Real,S<:Real}  

    E = cons.r*gₛ*VPD/P
    return E
end

#Calculate canopy transpiration (mol H₂O m⁻² leaf area s⁻¹)
function calc_Ec(E::Real,model::CCPHStruct)
    #E leaf transpiration (mol H₂O m⁻² leaf area s⁻¹)
    scaling_fac = calc_scaling_fac(model)
    Ec = E*scaling_fac
    return Ec
end

#Calculate soil-canopy conductance
function calc_K_cost(gₛ::S,
    H::T,
    hydPar::HydraulicsPar,
    env::EnvironmentStruct,
    cons::Constants;
    P_crit::Real=0.12) where {T<:Real,S<:Real}

    ψ₅₀,b,Kₓₗ₀,g,ρ_H2O,ψₛ = hydPar.ψ₅₀,hydPar.b,hydPar.Kₓₗ₀,cons.g,cons.ρ_H2O,env.ψₛ
    
    #Calculate transpiration  
    E = calc_E(gₛ,env.VPD,env.P;cons=cons)    

    #Soil water potential adjusted for gravitational pressure (MPa)
    ψₛ_g = ψₛ-H*ρ_H2O*g*10^-6  

    #Calculate canopy conductance    
    try      
        P = fixpoint(x->Ptarget(x,Kₓₗ₀,E,ψₛ_g,ψ₅₀,b),1.0)
        #f(x) = P_zero(x,Kₓₗ₀,E,ψₛ_g,ψ₅₀,b)
        #D(f) = x -> ForwardDiff.derivative(f,x)
        #P = Roots.find_zero((f, D(f)), 0.5, Roots.Newton();atol=0.0001)
        #P = bisection(f,0.0001,1.0)
        K_crit = Kₓₗ₀*P_crit
        Kₓₗ = Kₓₗ₀*P
        K_cost = (Kₓₗ-K_crit)/(Kₓₗ₀-K_crit)
        ψ_c = ψₛ_g-E/Kₓₗ 
        return K_cost, Kₓₗ, ψ_c     
    catch
        error("K_cost failed: gₛ=$(gₛ), E=$(E), Kₓₗ₀=$(Kₓₗ₀), ψₛ_g =$(ψₛ_g)")
    end      
end
function calc_K_cost(gₛ::T,model::CCPHStruct;P_crit::Real=0.12) where {T<:Real}
    
    #Calculate tree height
    H = model.treesize.H

    K_cost, Kₓₗ, ψ_c = calc_K_cost(gₛ,H,model.hydPar,model.env,model.cons;P_crit=P_crit) 

    return K_cost, Kₓₗ, ψ_c     
end

#Calculate the stomatal conductance (gₛ/E) such that Pfun(ψ_target,ψ₅₀,b)=Pval and E = Pint(ψ_target,ψₛ_g,ψ₅₀,b)
function calc_K_costᵢₙᵥ(Pval::Real,model::CCPHStruct)
    ψ₅₀,b,Kₓₗ₀,g,ρ_H2O,ψₛ = model.hydPar.ψ₅₀,model.hydPar.b,model.hydPar.Kₓₗ₀,model.cons.g,model.cons.ρ_H2O,model.env.ψₛ

    #Calc target ψ
    ψ_target = Pfunᵢₙᵥ(Pval,ψ₅₀,b)    
    #Calculate tree height
    H = model.treesize.H       
    #Soil water potential adjusted for gravitational pressure (MPa)
    ψₛ_g = ψₛ-H*ρ_H2O*g*10^-6     
    #Calculate target E
    E_target = Kₓₗ₀*Pint(ψ_target,ψₛ_g,ψ₅₀,b)
    #Calcualte target gₛ
    gₛ_target = E_target*model.env.P/(model.cons.r*model.env.VPD)

    return gₛ_target
end