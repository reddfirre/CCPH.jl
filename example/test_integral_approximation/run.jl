using CCPH, Plots
import BenchmarkTools

#Blahak 2010 approximation
function get_gamma_param(a::Real)
    p₁ = 9.4368392235E-3
    p₂ = -1.0782666481E-4
    p₃ = -5.8969657295E-6
    p₄ = 2.8939523781E-7
    p₅ = 1.0043326298E-1
    p₆ = 5.5637848465E-1

    q₁ = 1.1464706419E-1
    q₂ = 2.6963429121
    q₃ = -2.9647038257
    q₄ = 2.1080724954

    r₁ = 0.0
    r₂ = 1.1428716184
    r₃ = -6.6981186438E-3
    r₄ = 1.0480765092E-4

    s₁ = 1.0356711153
    s₂ = 2.3423452308
    s₃ = -3.6174503174E-1
    s₄ = -3.1376557650
    s₅ = 2.9092306039

    c₁ = 1+p₁*a+p₂*a^2+p₃*a^3+p₄*a^4+p₅*(exp(-p₆*a)-1)
    c₂ = q₁+q₂/a+q₃/a^2+q₄/a^3
    c₃ = r₁+r₂*a+r₃*a^2+r₄*a^3
    c₄ = s₁+s₂/a+s₃/a^2+s₄/a^3+s₅/a^4

    return c₁,c₂,c₃,c₄
end

function γₗ(a::Real,x::Real,c₁::Real,c₂::Real,c₃::Real,c₄::Real,Γₐ::Real)
    W=(1+tanh(c₂*(x-c₃)))/2
    cx = c₁*x
    ainv = 1/a
    aainv = ainv/(a+1)
    aaainv = aainv/(a+2)
    return exp(-x)*x^a*(ainv+cx*aainv+cx^2*aaainv)*(1-W)+Γₐ*W*(1-c₄^(-x))
end

Γᵤ(a::Real,x::Real,c₁::Real,c₂::Real,c₃::Real,c₄::Real,Γₐ::Real) = Γₐ-γₗ(a,x,c₁,c₂,c₃,c₄,Γₐ)

function CCPH.Pintlim(ψ::S,ψ₅₀::T,b::T,c₁::Real,c₂::Real,c₃::Real,c₄::Real,Γₐ::T) where {S<:Real,T<:Real} 
    Γ = Γᵤ((1+b)/b,log(2)*(ψ/ψ₅₀)^b,c₁,c₂,c₃,c₄,Γₐ)
    return CCPH.Pfun(ψ,ψ₅₀,b)*ψ+abs(ψ₅₀)/log(2)^(1/b)*Γ
end

function CCPH.Pint(ψ_c::S,ψ_cm::T,ψ₅₀::T,b::T,c₁::Real,c₂::Real,c₃::Real,c₄::Real,Γₐ::T) where {S<:Real,T<:Real} 
    return CCPH.Pintlim(ψ_cm,ψ₅₀,b,c₁,c₂,c₃,c₄,Γₐ)-CCPH.Pintlim(ψ_c,ψ₅₀,b,c₁,c₂,c₃,c₄,Γₐ)
end

function CCPH.Ptarget(x::W,Kₓₗ₀::T,E::S,ψ_cm::T,ψ₅₀::T,b::T,c₁::Real,c₂::Real,c₃::Real,c₄::Real,Γₐ::T) where {S<:Real,W<:Real,T<:Float64}
    return CCPH.Pint(ψ_cm-E/(Kₓₗ₀*x),ψ_cm,ψ₅₀,b,c₁,c₂,c₃,c₄,Γₐ)/(E/(Kₓₗ₀*x))
end

Pint_correct(ψ_c::S,ψ_cm::T,ψ₅₀::T,b::T) where {S<:Real,T<:Float64} = CCPH.Pintlim(ψ_cm,ψ₅₀,b)-CCPH.Pintlim(ψ_c,ψ₅₀,b)

function Ptarget_correct(x::W,Kₓₗ₀::T,E::S,ψ_cm::T,ψ₅₀::T,b::T) where {S<:Real,W<:Real,T<:Float64}
    return Pint_correct(ψ_cm-E/(Kₓₗ₀*x),ψ_cm,ψ₅₀,b)/(E/(Kₓₗ₀*x))
end

function Pint_Eller(ψ_c::Real,ψ_cm::Real,ψ₅₀::Real,b::Real)
    ψ = (ψ_c+ψ_cm)/2
    Δψ = (ψ_cm-ψ_c)
    return  CCPH.Pfun(ψ,ψ₅₀,b)*Δψ
end

#test integral approxiamtions
function run()
    ψₛ = CCPH.θₛ2ψₛ(0.06985236361)    
    H = 21.36 #m
    cons = Constants()
    ψₛ_g = ψₛ-H*cons.ρ_H2O*cons.g*10^-6 
    Kₓₗ₀ = 0.00075
    gₛ = 0.06    
    @show VPD = CCPH.VPDₜ(22.8,CCPH.SVPₜ(12.7))
    P = 1.0*10^5
    E = CCPH.calc_E(gₛ,VPD,P;cons=cons)        
    ψ₅₀ = -2.89
    b = 2.15
    @show ψ_crit = CCPH.Pfunᵢₙᵥ(0.12,ψ₅₀,b)
    @show E_crit = Kₓₗ₀*Pint_correct(ψ_crit,ψₛ_g,ψ₅₀,b)
    @show gₛ_crit = E_crit*P/(cons.r*VPD)
    @show ψ_min = CCPH.Pfunᵢₙᵥ(0.99,ψ₅₀,b)
    @show ψ_c = ψₛ_g-E/(Kₓₗ₀*0.8)
    @show a = (1+b)/b
    @show x = log(2)*(ψ_c/ψ₅₀)^b
    @show Γₐ = CCPH.SpecialFunctions.gamma(a)

    ψ_c_vec = range(ψ_min,stop=ψ_crit,length=200)

    #Correct integration values
    int_correct = Pint_correct.(ψ_c_vec,Ref(ψₛ_g),Ref(ψ₅₀),Ref(b))

    #Get parameters for the Blahak 2010 approximation
    c₁,c₂,c₃,c₄ = get_gamma_param(a)
    #Approximated integration using Blahak 2010 approximation
    int_Blahak = CCPH.Pint.(ψ_c_vec,Ref(ψₛ_g),Ref(ψ₅₀),Ref(b),Ref(c₁),Ref(c₂),Ref(c₃),Ref(c₄),Ref(Γₐ))

    #Approximated integration using Simpson's 1/3 rule
    int_simpson = CCPH.Pint.(ψ_c_vec,Ref(ψₛ_g),Ref(ψ₅₀),Ref(b))

    #Approximation used in Eller 2018
    int_Eller = Pint_Eller.(ψ_c_vec,Ref(ψₛ_g),Ref(ψ₅₀),Ref(b))

    pl1= plot(ψ_c_vec,int_correct,xlabel = "ψ_c", ylabel = "Integral",label="Exact")
    plot!(ψ_c_vec,int_Blahak,label="Blahak")
    plot!(ψ_c_vec,int_simpson,label="Simpson")
    plot!(ψ_c_vec,int_Eller,label="Eller")

    pl2 = plot(ψ_c_vec,(int_Blahak-int_correct)./int_correct*100,xlabel = "ψ_c", ylabel = "Relative Error %",label="Blahak")
    plot!(ψ_c_vec,(int_simpson-int_correct)./int_correct*100,label="Simpson")
    plot!(ψ_c_vec,(int_Eller-int_correct)./int_correct*100,label="Eller")

    plot(pl1,pl2,layout=(1,2))
end

run()