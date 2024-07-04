using CCPH, CSV, DataFrames, Statistics, StatsPlots
import BlackBoxOptim, ForwardDiff, BenchmarkTools

mutable struct OptVal
    gₛ₁::Real
    gₛ₂::Real
    Nₘ_f::Real
end

function CCPH.WeatherDataStruct(data::DataFrames.DataFrame,data_idx::Integer;lat::Real=64,Cₐ::Real=400.0/10.0,P::Real=1.0*10^5)
    d = data.Date[data_idx]
    data_day = CCPH.WeatherDataStruct(d,
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

function get_env_from_data(data::CCPH.WeatherDataStruct)
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

function get_env_from_data(data_vec::Vector{CCPH.WeatherDataStruct})
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

#Run optimization using global optimizer
function CCPHOpt_global(daylength::Real,
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

#Run optimization using global optimizer
function CCPHOpt_global(Nₘ_f::Real,daylength::Real,
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

function get_model_output_global(daylength::Real,model::CCPH.CCPHStruct,kinetic::CCPH.PhotoKineticRates,envfun::CCPH.EnvironmentFunStruct)
    gₛ₁_lim_hi,gₛ₂_lim_hi = CCPH.SDM2_get_gₛ_lim!(daylength,model,kinetic,envfun)     
    gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt = CCPHOpt_global(daylength,
    kinetic,
    envfun,
    model;
    gₛ₁_lim_hi=gₛ₁_lim_hi,
    gₛ₂_lim_hi=gₛ₂_lim_hi)
    optval = OptVal(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt)
    output = CCPH_run!(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt,daylength,kinetic,envfun,model)   
     
    return optval,output
end

function get_model_output_global(Nₘ_f_opt::Real,daylength::Real,model::CCPH.CCPHStruct,kinetic::CCPH.PhotoKineticRates,envfun::CCPH.EnvironmentFunStruct)
    gₛ₁_lim_hi,gₛ₂_lim_hi = CCPH.SDM2_get_gₛ_lim!(daylength,model,kinetic,envfun)     
    gₛ₁_opt,gₛ₂_opt = CCPHOpt_global(Nₘ_f_opt,
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

r_num(val::Real;digits::Integer=5) = round(val, digits=digits)

function comp_output(output::OptVal,output_alt::OptVal)
    println("Gradient: gₛ₁=$(r_num(output.gₛ₁)), gₛ₂=$(r_num(output.gₛ₂)), Nₘ_f=$(r_num(output.Nₘ_f))")
    println("Global: gₛ₁=$(r_num(output_alt.gₛ₁)), gₛ₂=$(r_num(output_alt.gₛ₂)), Nₘ_f=$(r_num(output_alt.Nₘ_f))")
    return nothing
end

function intitiate_model(data::DataFrames.DataFrame)
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

    data_week= CCPH.WeatherDataStruct.(Ref(data),1:7)

    env,envfun,daylength = get_env_from_data(data_week)    
    photo = PhotoPar(kinetic,env.Tₐ)
    model = CCPHStruct(cons,env,treepar,treesize,photo,hydPar)

    return model,envfun,daylength,kinetic,cons,treepar,treesize,hydPar
end

function intitiate_model(data::DataFrames.DataFrame,
    i::Integer,
    kinetic::CCPH.PhotoKineticRates,
    cons::CCPH.Constants,
    treepar::CCPH.TreePar,
    treesize::TreeSize,
    hydPar::HydraulicsPar)

    data_day = CCPH.WeatherDataStruct(data,i)
    env_day,envfun_day,daylength_day = get_env_from_data(data_day)    
    photo_day = PhotoPar(kinetic,env_day.Tₐ)
    model_day = CCPHStruct(cons,env_day,treepar,treesize,photo_day,hydPar)

    return model_day,envfun_day,daylength_day
end

function run_week(data::DataFrames.DataFrame)
    
    model,envfun,daylength,kinetic,cons,treepar,treesize,hydPar = intitiate_model(data)

    #Optimize weekly mean using BFGS
    optval_weekly,output_weekly = get_model_output(daylength,model,kinetic,envfun)
    
    for i = 1:7

        model_day,envfun_day,daylength_day = intitiate_model(data,i,kinetic,cons,treepar,treesize,hydPar)
        
        #Optimize daily using BFGS
        optval_day,output_day = get_model_output(optval_weekly.Nₘ_f,
        optval_weekly.gₛ₁,
        optval_weekly.gₛ₂,
        daylength_day,
        model_day,
        kinetic,
        envfun_day)           
    end
    
    return nothing
end

function run_week_global(data::DataFrames.DataFrame)

    model,envfun,daylength,kinetic,cons,treepar,treesize,hydPar = intitiate_model(data)

    #Optimize weekly mean using Global Opt
    optval_weekly_global,output_weekly_global = get_model_output_global(daylength,model,kinetic,envfun)    

    for i = 1:7

        model_day,envfun_day,daylength_day = intitiate_model(data,i,kinetic,cons,treepar,treesize,hydPar)
        
        #Optimize daily using Global Opt
        optval_day_global,output_day_global = get_model_output_global(optval_weekly_global.Nₘ_f,        
        daylength_day,
        model_day,
        kinetic,
        envfun_day)              
    end

    return nothing
end

function get_output_week(data::DataFrames.DataFrame)
    gₛ₁_output = zeros(8,2)
    gₛ₂_output = zeros(8,2)
    Nₘ_f_output = zeros(8,2)
    GPP_output = zeros(8,2)
    Ec_output = zeros(8,2)

    model,envfun,daylength,kinetic,cons,treepar,treesize,hydPar = intitiate_model(data)

    #Optimize weekly mean using BFGS
    optval_weekly,output_weekly = get_model_output(daylength,model,kinetic,envfun)

    #Optimize weekly mean using Global Opt
    optval_weekly_global,output_weekly_global = get_model_output_global(daylength,model,kinetic,envfun)
    
    gₛ₁_output[1,1] = optval_weekly.gₛ₁
    gₛ₁_output[1,2] = optval_weekly_global.gₛ₁

    gₛ₂_output[1,1] = optval_weekly.gₛ₂
    gₛ₂_output[1,2] = optval_weekly_global.gₛ₂

    Nₘ_f_output[1,1] = optval_weekly.Nₘ_f
    Nₘ_f_output[1,2] = optval_weekly_global.Nₘ_f

    GPP_output[1,1] = output_weekly.GPP
    GPP_output[1,2] = output_weekly_global.GPP

    Ec_output[1,1] = output_weekly.Ec
    Ec_output[1,2] = output_weekly_global.Ec

    for i = 1:7

        model_day,envfun_day,daylength_day = intitiate_model(data,i,kinetic,cons,treepar,treesize,hydPar)
        
        #Optimize daily using BFGS
        optval_day,output_day = get_model_output(optval_weekly.Nₘ_f,
        optval_weekly.gₛ₁,
        optval_weekly.gₛ₂,
        daylength_day,
        model_day,
        kinetic,
        envfun_day)    
        
        #Optimize daily using Global Opt
        optval_day_global,output_day_global = get_model_output_global(optval_weekly_global.Nₘ_f,        
        daylength_day,
        model_day,
        kinetic,
        envfun_day)

        gₛ₁_output[i+1,1] = optval_day.gₛ₁
        gₛ₁_output[i+1,2] = optval_day_global.gₛ₁
    
        gₛ₂_output[i+1,1] = optval_day.gₛ₂
        gₛ₂_output[i+1,2] = optval_day_global.gₛ₂
    
        Nₘ_f_output[i+1,1] = optval_day.Nₘ_f
        Nₘ_f_output[i+1,2] = optval_day_global.Nₘ_f
    
        GPP_output[i+1,1] = output_day.GPP
        GPP_output[i+1,2] = output_day_global.GPP
    
        Ec_output[i+1,1] = output_day.Ec
        Ec_output[i+1,2] = output_day_global.Ec
    end
    
    return gₛ₁_output,gₛ₂_output,Nₘ_f_output,GPP_output,Ec_output
end

function compare_output(data::DataFrames.DataFrame)
    
    gₛ₁_output,gₛ₂_output,Nₘ_f_output,GPP_output,Ec_output = get_output_week(data)

    ctg = repeat(["BFGS", "Global"], inner = 8)
    nam = repeat(["W", "D1", "D2", "D3", "D4", "D5", "D6", "D7"], outer = 2)
    
    ctg_Er = repeat(["Error"], inner = 8)
    nam_Er = ["W", "D1", "D2", "D3", "D4", "D5", "D6", "D7"]
    
    pl1 = groupedbar(nam,gₛ₁_output,group = ctg,ylabel="gₛ₁",
    bar_width = 0.67,lw = 0, framestyle = :box,yguidefontsize=8) 

    pl11 = groupedbar(nam_Er,transpose(100*abs.(gₛ₁_output[:,2]-gₛ₁_output[:,1])./gₛ₁_output[:,2]),group = ctg_Er,ylabel="Abs. Rel. Error % gₛ₁",
    bar_width = 0.67,lw = 0, framestyle = :box,legends=false,yguidefontsize=8)

    pl2 = groupedbar(nam,gₛ₂_output,group = ctg,ylabel="gₛ₂",
    bar_width = 0.67,lw = 0, framestyle = :box,yguidefontsize=8)

    pl22 = groupedbar(nam_Er,transpose(100*abs.(gₛ₂_output[:,2]-gₛ₂_output[:,1])./gₛ₂_output[:,2]),group = ctg_Er,ylabel="Abs. Rel. Error % gₛ₂",
    bar_width = 0.67,lw = 0, framestyle = :box,legends=false,yguidefontsize=8)

    pl3 = groupedbar(nam,Nₘ_f_output,group = ctg,ylabel="Nₘ_f",
    bar_width = 0.67,lw = 0, framestyle = :box,yguidefontsize=8)

    pl33 = groupedbar(nam_Er,transpose(100*abs.(Nₘ_f_output[:,2]-Nₘ_f_output[:,1])./Nₘ_f_output[:,2]),group = ctg_Er,ylabel="Abs. Rel. Error % Nₘ_f",
    bar_width = 0.67,lw = 0, framestyle = :box,legends=false,yguidefontsize=8)   

    plot(pl1,pl11,pl2,pl22,pl3,pl33,layout=(3,2))
end

#Test global optimizer vs gradient-based optimizer
function run_time_test()   
    #Check run-time
    for i = 1:4
        println("--Test week $(i)--")
        #Load data
        data = CSV.read("Test_Week_$(i)_Data.csv", DataFrame)
        
        BenchmarkTools.@btime run_week($data)    
        BenchmarkTools.@btime run_week_global($data)  
    end
end  

function opt_val_test()
    week_idx = 4
    #Load data    
    data = CSV.read("Test_Week_$(week_idx)_Data.csv", DataFrame)
    #Display optimal values (BFGS vs Global Opt)
    compare_output(data::DataFrames.DataFrame)
end

run_time_test()
opt_val_test()