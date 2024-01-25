using CCPH, CSV, DataFrames, Statistics
import BlackBoxOptim, ForwardDiff, BenchmarkTools

mutable struct OptVal
    gₛ₁::Real
    gₛ₂::Real
    Nₘ_f::Real
end

mutable struct DataStruct
    date::CCPH.Dates.Date #Date yyyy-mm-dd
    lat::Real #Latitude in degrees
    Tmean::Real #°C
    Tmin::Real #°C
    Tmax::Real #°C
    VP::Real #Pa
    Radₜₒ::Real #J m⁻²
    θₛ::Real #-
    Cₐ::Real #Pa
    P::Real #Pa
end
function DataStruct(data::DataFrames.DataFrame,data_idx::Integer;lat::Real=64,Cₐ::Real=400.0/10.0,P::Real=1.0*10^5)
    #d = CCPH.Dates.Date(data.Date[data_idx], "dd/mm/yyyy")
    d = data.Date[data_idx]
    data_day = DataStruct(d,
    lat,
    data.airTmean[data_idx],
    data.airTmin[data_idx],
    data.airTmax[data_idx]
    ,data.VP[data_idx],
    data.Radiation[data_idx]*10^6,
    data.SWC[data_idx],
    Cₐ,
    P)
    return data_day
end

function get_env_from_data(data::DataStruct)
    day_nr = CCPH.Dates.dayofyear(data.date)
    daylength = CCPH.daylighthour(data.lat*pi/180,day_nr)*3600 #Seconds
    I₀ₜₒₜ = data.Radₜₒ*2.3*10^-6 #mol m⁻²    
    VP_Tmin = CCPH.SVPₜ(data.Tmin) #Pa
    VP = min(VP_Tmin,data.VP)

    ψₛ = CCPH.θₛ2ψₛ(data.θₛ) #MPa
    Cₐ = data.Cₐ #Pa
    P = data.P #Pa
    Tₐ = data.Tmean #°C
    VPD =  CCPH.VPDₜ(Tₐ,VP) #Pa
    I₀ = I₀ₜₒₜ/daylength #mol m⁻² s⁻¹

    ψₛ_fun(t) = ψₛ
    Cₐ_fun(t) = Cₐ
    P_fun(t) = P
    Tₐ_fun(t) = CCPH.Tₐ_fun(t,data.Tmin,data.Tmax,daylength)
    VPD_fun(t) = CCPH.VPDₜ(Tₐ_fun(t),VP)
    I₀_fun(t) = CCPH.I₀_fun(t,I₀ₜₒₜ,daylength)

    envfun = EnvironmentFunStruct(I₀_fun,Cₐ_fun,P_fun,Tₐ_fun,VPD_fun,ψₛ_fun)
    env = EnvironmentStruct(I₀,Cₐ,P,Tₐ,VPD,ψₛ)

    return env,envfun,daylength
end

function get_env_from_data(data_vec::Vector{DataStruct})
    day_nrs = [CCPH.Dates.dayofyear(data.date) for data in data_vec]
    lats = [data.lat for data in data_vec]
    daylengths = CCPH.daylighthour.(lats*pi/180,day_nrs)*3600 #Seconds
    daylength = mean(daylengths)

    I₀ₜₒₜ = mean([data.Radₜₒ for data in data_vec])*2.3*10^-6 #mol m⁻²  
    Tmin_data = mean([data.Tmin for data in data_vec]) #°C
    Tmax_data = mean([data.Tmax for data in data_vec]) #°C

    VP_Tmin = CCPH.SVPₜ(Tmin_data) #Pa  
    VP_data = mean([data.VP for data in data_vec])
    VP = min(VP_Tmin,VP_data)

    ψₛ = mean([CCPH.θₛ2ψₛ(data.θₛ) for data in data_vec]) #MPa
    Cₐ = mean([data.Cₐ for data in data_vec]) #Pa
    P =  mean([data.P for data in data_vec]) #Pa
    Tₐ =  mean([data.Tmean for data in data_vec]) #°C

    VPD =  CCPH.VPDₜ(Tₐ,VP) #Pa
    I₀ = I₀ₜₒₜ/daylength #mol m⁻² s⁻¹

    ψₛ_fun(t) = ψₛ
    Cₐ_fun(t) = Cₐ
    P_fun(t) = P
    Tₐ_fun(t) = CCPH.Tₐ_fun(t,Tmin_data,Tmax_data,daylength)
    VPD_fun(t) = CCPH.VPDₜ(Tₐ_fun(t),VP)
    I₀_fun(t) = CCPH.I₀_fun(t,I₀ₜₒₜ,daylength)

    envfun = EnvironmentFunStruct(I₀_fun,Cₐ_fun,P_fun,Tₐ_fun,VPD_fun,ψₛ_fun)
    env = EnvironmentStruct(I₀,Cₐ,P,Tₐ,VPD,ψₛ)

    return env,envfun,daylength
end

function CCPHOpt_alt(daylength::Real,
    photo_kinetic::CCPH.PhotoKineticRates,
    envfun::CCPH.EnvironmentFunStruct,
    model::CCPH.CCPHStruct;
    gₛ₁_lim_lo::Real=0.001,gₛ₁_lim_hi::Real=0.5,
    gₛ₂_lim_lo::Real=0.001,gₛ₂_lim_hi::Real=0.5,
    Nₘ_f_lim_lo::Real=0.001,Nₘ_f_lim_hi::Real=0.05,
    P_crit::Real=0.12)    
    
    res = BlackBoxOptim.bboptimize(x->CCPH.Objective_fun(x,daylength,photo_kinetic,envfun,model;P_crit=P_crit); SearchRange =[(gₛ₁_lim_lo,gₛ₁_lim_hi),(gₛ₂_lim_lo,gₛ₂_lim_hi),(Nₘ_f_lim_lo,Nₘ_f_lim_hi)],TraceMode=:silent)
    gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt = BlackBoxOptim.best_candidate(res)    
    
    return gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt
end

function CCPHOpt_alt(Nₘ_f::Real,daylength::Real,
    photo_kinetic::CCPH.PhotoKineticRates,
    envfun::CCPH.EnvironmentFunStruct,
    model::CCPH.CCPHStruct;
    gₛ₁_lim_lo::Real=0.001,gₛ₁_lim_hi::Real=0.5,
    gₛ₂_lim_lo::Real=0.001,gₛ₂_lim_hi::Real=0.5,
    P_crit::Real=0.12)    
    
    res = BlackBoxOptim.bboptimize(x->CCPH.Objective_fun([x[1],x[2],Nₘ_f],daylength,photo_kinetic,envfun,model;P_crit=P_crit); SearchRange =[(gₛ₁_lim_lo,gₛ₁_lim_hi),(gₛ₂_lim_lo,gₛ₂_lim_hi)],TraceMode=:silent)
    gₛ₁_opt,gₛ₂_opt = BlackBoxOptim.best_candidate(res)    
    
    return gₛ₁_opt,gₛ₂_opt
end

function get_model_output(daylength::Real,model::CCPH.CCPHStruct,kinetic::CCPH.PhotoKineticRates,envfun::CCPH.EnvironmentFunStruct)
    gₛ₁_lim_hi,gₛ₂_lim_hi = CCPH.SDM2_get_gₛ_lim!(daylength,model,kinetic,envfun)     
    gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt = CCPHOpt(daylength,kinetic,envfun,model;gₛ₁_lim_hi=gₛ₁_lim_hi,gₛ₂_lim_hi=gₛ₂_lim_hi)
    optval = OptVal(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt)
    output = CCPH_run!(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt,daylength,kinetic,envfun,model)
    return optval,output
end

function get_model_output(Nₘ_f_opt::Real,
    gₛ₁_guess::Real,
    gₛ₂_guess::Real,
    daylength::Real,
    model::CCPH.CCPHStruct,
    kinetic::CCPH.PhotoKineticRates,
    envfun::CCPH.EnvironmentFunStruct)

    gₛ₁_lim_hi,gₛ₂_lim_hi = CCPH.SDM2_get_gₛ_lim!(daylength,model,kinetic,envfun)     
    gₛ₁_opt,gₛ₂_opt = CCPHOpt(Nₘ_f_opt,
    daylength,
    kinetic,
    envfun,
    model;
    gₛ₁_lim_hi=gₛ₁_lim_hi,
    gₛ₂_lim_hi=gₛ₂_lim_hi,
    gₛ₁_guess=gₛ₁_guess,
    gₛ₂_guess=gₛ₂_guess)    
    optval = OptVal(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt)
    output = CCPH_run!(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt,daylength,kinetic,envfun,model)
    return optval,output
end

function get_model_output_alt(daylength::Real,model::CCPH.CCPHStruct,kinetic::CCPH.PhotoKineticRates,envfun::CCPH.EnvironmentFunStruct)
    gₛ₁_lim_hi,gₛ₂_lim_hi = CCPH.SDM2_get_gₛ_lim!(daylength,model,kinetic,envfun)     
    gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt = CCPHOpt_alt(daylength,
    kinetic,
    envfun,
    model;
    gₛ₁_lim_hi=gₛ₁_lim_hi,
    gₛ₂_lim_hi=gₛ₂_lim_hi)
    optval = OptVal(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt)
    output = CCPH_run!(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt,daylength,kinetic,envfun,model)   
     
    return optval,output
end

function get_model_output_alt(Nₘ_f_opt::Real,daylength::Real,model::CCPH.CCPHStruct,kinetic::CCPH.PhotoKineticRates,envfun::CCPH.EnvironmentFunStruct)
    gₛ₁_lim_hi,gₛ₂_lim_hi = CCPH.SDM2_get_gₛ_lim!(daylength,model,kinetic,envfun)     
    gₛ₁_opt,gₛ₂_opt = CCPHOpt_alt(Nₘ_f_opt,
    daylength,
    kinetic,
    envfun,
    model;
    gₛ₁_lim_hi=gₛ₁_lim_hi,
    gₛ₂_lim_hi=gₛ₂_lim_hi)
    optval = OptVal(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt)
    output = CCPH_run!(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt,daylength,kinetic,envfun,model)   
     
    return optval,output
end

r_num(val::Real;digits::Integer=4) = round(val, digits=digits)

function comp_output(output::OptVal,output_alt::OptVal)
    println("BFGS: gₛ₁=$(r_num(output.gₛ₁)), gₛ₂=$(r_num(output.gₛ₂)), Nₘ_f=$(r_num(output.Nₘ_f))")
    println("Global: gₛ₁=$(r_num(output_alt.gₛ₁)), gₛ₂=$(r_num(output_alt.gₛ₂)), Nₘ_f=$(r_num(output_alt.Nₘ_f))")
    return nothing
end

function run_week(data::DataFrames.DataFrame)
    #Stand size data
    H = 19.07 #m
    N = 850.0/10000 # #tree/m²
    LAI = 2.45 #-

    #Create structs
    cons = Constants()   
    treepar = TreePar(Nₛ=0.013,α_max=0.14,a_Jmax=0.011,b_Jmax=0.0)
    treesize = TreeSize(H,LAI,N)
    hydPar = HydraulicsPar(Kₓₗ₀=0.00054)
    kinetic = PhotoKineticRates()

    data_week= DataStruct.(Ref(data),1:7)

    env,envfun,daylength = get_env_from_data(data_week)    
    photo = PhotoPar(kinetic,env.Tₐ)
    model = CCPHStruct(cons,env,treepar,treesize,photo,hydPar)

    #Optimize weekly mean using BFGS
    optval_weekly,output_weekly = get_model_output(daylength,model,kinetic,envfun)
    #Optimize weekly mean using Global Opt
    #optval_weekly_alt,output_weekly_alt = get_model_output(daylength,model,kinetic,envfun)

    println("Weekly")
    #comp_output(optval_weekly,optval_weekly_alt)

    for i = 1:7
        println("Day nr $(i)")
        #Optimize daily

        data_day = DataStruct(data,i)
        env_day,envfun_day,daylength_day = get_env_from_data(data_day)    
        photo_day = PhotoPar(kinetic,env_day.Tₐ)
        model_day = CCPHStruct(cons,env_day,treepar,treesize,photo_day,hydPar)
        
        #Optimize daily using BFGS
        optval_day,output_day = get_model_output(optval_weekly.Nₘ_f,
        optval_weekly.gₛ₁,
        optval_weekly.gₛ₂,
        daylength_day,
        model_day,
        kinetic,
        envfun_day) 
        
        #=
        #Optimize daily using Global Opt
        optval_day_alt,output_day_alt = get_model_output_alt(optval_weekly_alt.Nₘ_f,        
        daylength_day,
        model_day,
        kinetic,
        envfun_day)
        =#  
        
        #comp_output(optval_day,optval_day_alt)
    end
    
    return nothing
end

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

function Pint_simps(ψ_c::S,ψ_cm::T,ψ₅₀::T,b::T)  where {S<:Real,T<:Real}
    ψ = (ψ_c+ψ_cm)/2
    Δψ = (ψ_cm-ψ_c)
    return Δψ*(CCPH.Pfun(ψ_c,ψ₅₀,b)+4*CCPH.Pfun(ψ,ψ₅₀,b)+CCPH.Pfun(ψ_cm,ψ₅₀,b))/6
end

function Ptarget_simps(x::W,Kₓₗ₀::T,E::S,ψ_cm::T,ψ₅₀::T,b::T) where {S<:Real,W<:Real,T<:Float64}
    return Pint_simps(ψ_cm-E/(Kₓₗ₀*x),ψ_cm,ψ₅₀,b)/(E/(Kₓₗ₀*x))
end

#Test global optimizer vs gradient-based optimizer
function run()
    #=
    ψₛ = CCPH.θₛ2ψₛ(0.2573929584)
    H = 19.07 #m
    cons = Constants()
    ψₛ_g = ψₛ-H*cons.ρ_H2O*cons.g*10^-6 
    Kₓₗ₀ = 0.00054
    gₛ = 0.06
    @show VPD = CCPH.VPDₜ(5.2,CCPH.SVPₜ(0.0))
    P = 1.0*10^5
    E = CCPH.calc_E(gₛ,VPD,P;cons=cons)        
    ψ₅₀ = -2.89
    b = 2.15
    @show ψ_c = ψₛ_g-E/(Kₓₗ₀*0.8)
    @show a = (1+b)/b
    @show x = log(2)*(ψ_c/ψ₅₀)^b
    @show Γₐ = CCPH.SpecialFunctions.gamma(a)
    
    p,q = CCPH.SpecialFunctions.gamma_inc(a,x)
    @show p*Γₐ

    c₁,c₂,c₃,c₄ = get_gamma_param(a)
    γₗ(a,x,c₁,c₂,c₃,c₄,Γₐ)

    @show CCPH.SpecialFunctions.gamma(a,x)
    Γᵤ(a,x,c₁,c₂,c₃,c₄,Γₐ)

    f() = CCPH.fixpoint(x->CCPH.Ptarget(x,Kₓₗ₀,E,ψₛ_g,ψ₅₀,b),1.0)
    f_alt() = CCPH.fixpoint(x->CCPH.Ptarget(x,Kₓₗ₀,E,ψₛ_g,ψ₅₀,b,Γₐ),1.0)
    f_alt2() = CCPH.fixpoint(x->CCPH.Ptarget(x,Kₓₗ₀,E,ψₛ_g,ψ₅₀,b,c₁,c₂,c₃,c₄,Γₐ),1.0)
    f_alt3() = CCPH.fixpoint(x->Ptarget_simps(x,Kₓₗ₀,E,ψₛ_g,ψ₅₀,b),1.0)

    BenchmarkTools.@btime $f()
    BenchmarkTools.@btime $f_alt()
    BenchmarkTools.@btime $f_alt2()
    BenchmarkTools.@btime $f_alt3()

    @show f()
    @show f_alt3()
    =#
    
    for i = 1:4
        println("--$(i)--")
        #Load data
        data = CSV.read("Test_Week_$(i)_Data.csv", DataFrame)
        
        @time run_week(data)    
    end    
end

run()