using CCPH, CSV, DataFrames, Statistics
import BlackBoxOptim, ForwardDiff

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
    d = CCPH.Dates.Date(data.Date[data_idx], "dd/mm/yyyy")
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

function CCPH.CCPHOpt(Nₘ_f::Real,daylength::Real,
    photo_kinetic::CCPH.PhotoKineticRates,
    envfun::CCPH.EnvironmentFunStruct,
    model::CCPH.CCPHStruct;
    gₛ₁_guess::Real=0.02,gₛ₁_lim_lo::Real=0.001,gₛ₁_lim_hi::Real=0.5,
    gₛ₂_guess::Real=0.02,gₛ₂_lim_lo::Real=0.001,gₛ₂_lim_hi::Real=0.5,
    P_crit::Real=0.12)
    
    x0 = [gₛ₁_guess, gₛ₂_guess]    
  
    lower = [gₛ₁_lim_lo, gₛ₂_lim_lo]   
    upper = [gₛ₁_lim_hi, gₛ₂_lim_hi] 
    
    res = BlackBoxOptim.bboptimize(x->CCPH.Objective_fun([x[1],x[2],Nₘ_f],daylength,photo_kinetic,envfun,model;P_crit=P_crit); SearchRange =[(gₛ₁_lim_lo,gₛ₁_lim_hi),(gₛ₂_lim_lo,gₛ₂_lim_hi)])
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

function get_model_output(Nₘ_f_opt::Real,daylength::Real,model::CCPH.CCPHStruct,kinetic::CCPH.PhotoKineticRates,envfun::CCPH.EnvironmentFunStruct)
    gₛ₁_lim_hi,gₛ₂_lim_hi = CCPH.SDM2_get_gₛ_lim!(daylength,model,kinetic,envfun)     
    gₛ₁_opt,gₛ₂_opt = CCPHOpt(Nₘ_f_opt,daylength,kinetic,envfun,model;gₛ₁_lim_hi=gₛ₁_lim_hi,gₛ₂_lim_hi=gₛ₂_lim_hi)
    optval = OptVal(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt)
    output = CCPH_run!(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt,daylength,kinetic,envfun,model)
    return optval,output
end
function g_fun!(storage, x,f)

        df = ForwardDiff.gradient(f, x)

        storage[1] = df[1]
        storage[2] = df[2]        
end

function CCPHOpt_alt(Nₘ_f::Real,daylength::Real,
    photo_kinetic::CCPH.PhotoKineticRates,
    envfun::CCPH.EnvironmentFunStruct,
    model::CCPH.CCPHStruct;
    gₛ₁_guess::Real=0.02,gₛ₁_lim_lo::Real=0.001,gₛ₁_lim_hi::Real=0.5,
    gₛ₂_guess::Real=0.02,gₛ₂_lim_lo::Real=0.001,gₛ₂_lim_hi::Real=0.5,
    P_crit::Real=0.12)

    x0 = [gₛ₁_guess, gₛ₂_guess]    
  
    lower = [gₛ₁_lim_lo, gₛ₂_lim_lo]   
    upper = [gₛ₁_lim_hi, gₛ₂_lim_hi] 

    
    
    #f(x) = CCPH.Objective_fun([x[1],x[2],Nₘ_f],daylength,photo_kinetic,envfun,model;P_crit=P_crit)
    #g!(storage, x) = g_fun!(storage, x,f)    

    #df = CCPH.Optim.OnceDifferentiable(f,g!,x0)   
    
    #inner_optimizer = CCPH.Optim.BFGS()
    #opt_trait = CCPH.Optim.optimize(df, lower, upper, x0, CCPH.Optim.Fminbox(inner_optimizer))
    
  
    CCPH.Optim.converged(opt_trait)||error("No optimal traits could be found") 
      
    
    gₛ₁_opt = opt_trait.minimizer[1]
    gₛ₂_opt = opt_trait.minimizer[2]

    return gₛ₁_opt,gₛ₂_opt
end
function get_model_output(Nₘ_f_opt::Real,gₛ₁_guess::Real,gₛ₂_guess::Real,daylength::Real,model::CCPH.CCPHStruct,kinetic::CCPH.PhotoKineticRates,envfun::CCPH.EnvironmentFunStruct)
    gₛ₁_lim_hi,gₛ₂_lim_hi = CCPH.SDM2_get_gₛ_lim!(daylength,model,kinetic,envfun)     
    gₛ₁_opt,gₛ₂_opt = CCPHOpt_alt(Nₘ_f_opt,
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

function run()
    #Stand size data
    H = 19.07 #m
    N = 850.0/10000 # #tree/m²
    LAI = 2.45 #-

    #Load data
    data = CSV.read("Test_Week_Data.csv", DataFrame)
    println(data)

    #data_day = DataStruct(data,1)   

    #env,envfun,daylength = get_env_from_data(data_day)

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

    #Optimize weekly mean
    optval_weekly,output_weekly = get_model_output(daylength,model,kinetic,envfun)
    @show optval_weekly

    #Optimize daily
    data_day = DataStruct(data,1)
    env_day,envfun_day,daylength_day = get_env_from_data(data_day)    
    photo_day = PhotoPar(kinetic,env_day.Tₐ)
    model_day = CCPHStruct(cons,env_day,treepar,treesize,photo_day,hydPar)
    optval_day,output_day = get_model_output(optval_weekly.Nₘ_f,daylength_day,model_day,kinetic,envfun_day)
    optval_day_altopt,output_day_altopt = get_model_output(optval_weekly.Nₘ_f,optval_weekly.gₛ₁,optval_weekly.gₛ₂,daylength_day,model_day,kinetic,envfun_day)
    @show optval_day

    output_day_alt = CCPH_run!(optval_weekly.gₛ₁,optval_weekly.gₛ₂,optval_weekly.Nₘ_f,daylength_day,kinetic,envfun_day,model_day)

    @show output_day
    @show output_day_alt    
    @show output_day_altopt
end

run()