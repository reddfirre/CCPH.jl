using CCPH, CSV, DataFrames, Plots, Statistics

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
    ,data.VP[data_idx]*100,
    data.Radiation[data_idx]*10^6,
    data.SWC[data_idx]/100,
    Cₐ,
    P)
    return data_day
end

#The delayed temperature, Sₜ, is modelled using a first order delay dynamics model
Sₜ_fun(Tₜ::Real,Sₜ₋₁::Real,τ::Real) = (1-1/τ)*Sₜ₋₁+Tₜ/τ
function Sₜ_fun(data::Vector{CCPH.WeatherDataStruct},τ::Real)
    n_data = length(data)
    S₁ = data[1].Tmean
    S_vec = [S₁]
    for i in 2:n_data
        Tₜ = data[i].Tmean   
        Sₜ₋₁ = S_vec[i-1]
        push!(S_vec,Sₜ_fun(Tₜ,Sₜ₋₁,τ))
    end
    return S_vec
end

#Photosynthetic temperature acclimation factor Xₜ (Mäkelä et al., 2004; Mäkelä et al., 2008).
function Xₜ_fun(Sₜ::Real,Smin::Real,ΔS::Real)
    if Sₜ≤Smin
        return 0
    elseif Smin<Sₜ<Smin+ΔS
        return (Sₜ-Smin)/ΔS
    elseif Sₜ≥Smin+ΔS
        return 1
    else
        error("Sₜ")
    end
end
function Xₜ_fun(raw_data::Vector{CCPH.WeatherDataStruct};Smin::Real = -4.0,ΔS::Real = 16.0,τ::Real =  7.0)
    #Calcualte photosynthetic temperature acclimation factor Xₜ (Mäkelä et al., 2004; Mäkelä et al., 2008). 
    #Smin = -4.0 #minium tempreture for Photosynthesis    
    #τ =  7.0 #Days
    #Smax = 16.0  
    
    S_vec = Sₜ_fun(raw_data,τ)
    X_vec = Xₜ_fun.(S_vec,Ref(Smin),Ref(ΔS))
    return X_vec
end

meanVPD(data::CCPH.WeatherDataStruct) = CCPH.VPDₜ(data.Tmean,data.VP)
get_daylength(data::CCPH.WeatherDataStruct) = CCPH.daylighthour(data.lat*pi/180,CCPH.Dates.dayofyear(data.date))*3600


function calc_Ec_data(VPD_z::Real,Sₑ::Real;
    E_cm_ref::Real=1.812,s_VPD_z::Real=3.121,s_Sₑ::Real=18.342)
    #Tor-Ngern et al. 2017 model for estimating canopy transpiration E_c (mm/d) 
    #Default parameter values taken are form the original paper
    #for estimating E_c for Rosinedal Scots pine stand (both fertilized and control)
    #---Input---
    #VPD_z day-length normalized vapor pressure deficit (VPD*Day length/24, Pa) 
    #Sₑ effective saturation (-)
    E_cm = E_cm_ref*(1-exp(-s_VPD_z*VPD_z*10^-3))
    Ec = E_cm*(1-exp(-s_Sₑ*Sₑ))

    return Ec
end
function calc_Ec_data(data::CCPH.WeatherDataStruct;
    E_cm_ref::Real=1.812,
    s_VPD_z::Real=3.121,
    s_Sₑ::Real=18.342)

    VPD = meanVPD(data)
    daylength = get_daylength(data)    
    θₛ = data.θₛ       
    VPD_z = VPD*daylength/(24*3600)    
    Sₑ = CCPH.calcSₑ(θₛ)

    return calc_Ec_data(VPD_z,Sₑ;E_cm_ref=E_cm_ref,s_VPD_z=s_VPD_z,s_Sₑ=s_Sₑ)
end
function calc_Ec_data(data::Vector{CCPH.WeatherDataStruct};
    E_cm_ref::Real=1.812,
    s_VPD_z::Real=3.121,
    s_Sₑ::Real=18.342)

    return calc_Ec_data.(data;E_cm_ref=E_cm_ref,s_VPD_z=s_VPD_z,s_Sₑ=s_Sₑ)
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

function get_model_output(daylength::Real,model::CCPH.CCPHStruct,kinetic::CCPH.PhotoKineticRates,envfun::CCPH.EnvironmentFunStruct)
    gₛ₁_lim_hi,gₛ₂_lim_hi = CCPH.SDM2_get_gₛ_lim!(daylength,model,kinetic,envfun)     
    gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt = CCPH.CCPHOpt(daylength,kinetic,envfun,model;gₛ₁_lim_hi=gₛ₁_lim_hi,gₛ₂_lim_hi=gₛ₂_lim_hi)
    optval = OptVal(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt)
    output = CCPH.CCPH_run!(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt,daylength,kinetic,envfun,model)
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
    gₛ₁_opt,gₛ₂_opt = CCPH.CCPHOpt(Nₘ_f_opt,
    daylength,
    kinetic,
    envfun,
    model;
    gₛ₁_lim_hi=gₛ₁_lim_hi,
    gₛ₂_lim_hi=gₛ₂_lim_hi,
    gₛ₁_guess=min(gₛ₁_guess,0.9*gₛ₁_lim_hi),
    gₛ₂_guess=min(gₛ₂_guess,0.9*gₛ₂_lim_hi))    
    optval = OptVal(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt)
    output = CCPH.CCPH_run!(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt,daylength,kinetic,envfun,model)
    return optval,output
end

function intitiate_model(data_week::CCPH.WeatherDataStruct,Xₜ::Real)
    #Stand size data
    H = 19.07 #m
    N = 850.0/10000 # #tree/m²
    LAI = 2.45 #-

    #Create structs
    cons = Constants()   
    treepar = TreePar(Xₜ=Xₜ,Nₛ=0.013,α_max=0.14,a_Jmax=0.011,b_Jmax=0.0)
    treesize = TreeSize(H,LAI,N)
    hydPar = HydraulicsPar(Kₓₗ₀=0.00054)
    kinetic = PhotoKineticRates()    

    env,envfun,daylength = get_env_from_data(data_week)    
    photo = PhotoPar(kinetic,env.Tₐ)
    model = CCPHStruct(cons,env,treepar,treesize,photo,hydPar)

    return model,envfun,daylength,kinetic,cons,treepar,treesize,hydPar
end

function intitiate_model(data_week::Vector{CCPH.WeatherDataStruct},Xₜ::Vector{T}) where {T<:Real}
    #Stand size data
    H = 19.07 #m
    N = 850.0/10000 # #tree/m²
    LAI = 2.45 #-

    #Create structs
    cons = Constants()   
    treepar = TreePar(Xₜ=mean(Xₜ),Nₛ=0.013,α_max=0.14,a_Jmax=0.011,b_Jmax=0.0)
    treesize = TreeSize(H,LAI,N)
    hydPar = HydraulicsPar(Kₓₗ₀=0.00054)
    kinetic = PhotoKineticRates()    

    env,envfun,daylength = get_env_from_data(data_week)    
    photo = PhotoPar(kinetic,env.Tₐ)
    model = CCPHStruct(cons,env,treepar,treesize,photo,hydPar)

    return model,envfun,daylength,kinetic,cons,treepar,treesize,hydPar
end


function intitiate_model(data_week::Vector{CCPH.WeatherDataStruct},
    Xₜ_week::Vector{T},
    i::Integer,
    kinetic::CCPH.PhotoKineticRates,
    cons::CCPH.Constants,
    treepar::CCPH.TreePar,
    treesize::TreeSize,
    hydPar::HydraulicsPar) where {T<:Real}

    data_day = data_week[i] 
    Xₜ = Xₜ_week[i] 
    env_day,envfun_day,daylength_day = get_env_from_data(data_day)    
    photo_day = PhotoPar(kinetic,env_day.Tₐ)
    treepar.Xₜ = Xₜ
    model_day = CCPHStruct(cons,env_day,treepar,treesize,photo_day,hydPar)

    return model_day,envfun_day,daylength_day
end


function run_week(data::Vector{CCPH.WeatherDataStruct},Xₜ::Vector{T}) where {T<:Real}
    
    model,envfun,daylength,kinetic,cons,treepar,treesize,hydPar = intitiate_model(data,Xₜ)

    #Optimize weekly mean using BFGS
    optval_weekly,output_weekly = get_model_output(daylength,model,kinetic,envfun)
    
    optval_day_vec = Vector{OptVal}(undef,7)
    output_day_vec = Vector{CCPH.CCPHOutput}(undef,7)
    
    for i = 1:7

        model_day,envfun_day,daylength_day = intitiate_model(data,Xₜ,i,kinetic,cons,treepar,treesize,hydPar)
        
        #Optimize daily using BFGS
        optval_day,output_day = get_model_output(optval_weekly.Nₘ_f,
        optval_weekly.gₛ₁,
        optval_weekly.gₛ₂,
        daylength_day,
        model_day,
        kinetic,
        envfun_day)        
        
        optval_day_vec[i] = optval_day
        output_day_vec[i] = output_day
    end
    
    return optval_weekly,output_weekly,optval_day_vec,output_day_vec
end

function Load_RO_weather_data(year::Integer;stand_type::Symbol=:Fertilized)
    SWC_data = CSV.read("./Daily_SWC_data_$(year).csv", DataFrames.DataFrame)
    Weather_data = CSV.read("./Weather_data_$(year).csv", DataFrames.DataFrame)

    if stand_type==:Fertilized
        Weather_data[!,:SWC] = SWC_data.SWRos2         
    elseif stand_type==:Control
        Weather_data[!,:SWC] = SWC_data.SWRos3
    else           
        error("Wrong input. Either :Fertilized or :Control") 
    end

    n_rows =  DataFrames.nrow(Weather_data)
    raw_data = [CCPH.WeatherDataStruct(Weather_data,i) for i in 1:n_rows]

    return raw_data
end 
function get_GPP_data(year::Integer;stand_type::Symbol=:Fertilized)
    GPP_data = CSV.read("./GPP_data_$(year).csv", DataFrames.DataFrame)
 
     if stand_type==:Fertilized
         return GPP_data.GPP_RO2       
     elseif stand_type==:Control
         return GPP_data.GPP_RO3
     else           
         error("Wrong input. Either :Fertilized or :Control") 
     end
 end

 function old_obj_fun(x,model::CCPH.CCPHStruct)
    gₛ,Nₘ_f = x
    output = CCPH.CCPH_inst_run(gₛ,Nₘ_f,model)
    return -output.Gain
end

function run_week_old(data::Vector{CCPH.WeatherDataStruct},Xₜ::Vector{T}) where {T<:Real}    
    model,envfun,daylength,kinetic,cons,treepar,treesize,hydPar = intitiate_model(data,Xₜ)

    gₛ_lim_hi = CCPH.calc_K_costᵢₙᵥ(0.12,model)
    x0 = [0.02, 0.012]
    lower = [0.001, 0.001]   
    upper = [gₛ_lim_hi, 0.05] 

    df = CCPH.Optim.OnceDifferentiable(x->old_obj_fun(x,model),x0)   

    inner_optimizer = CCPH.Optim.BFGS(linesearch = CCPH.Optim.LineSearches.BackTracking())
    opt_trait = CCPH.Optim.optimize(df, lower, upper, x0, CCPH.Optim.Fminbox(inner_optimizer))

    CCPH.Optim.converged(opt_trait)||error("No optimal traits could be found") 

    gₛ_opt = opt_trait.minimizer[1]
    Nₘ_f_opt = opt_trait.minimizer[2]

    #Output from old model
    output = CCPH.CCPH_inst_run(gₛ_opt,Nₘ_f_opt,model)
    Ec = output.Ec*daylength
    Ec *= 1000*model.cons.M_H2O/model.cons.ρ_H2O
    
    GPP = output.GPP*daylength
    GPP *= model.cons.M_C*1000
    
    return GPP,Ec
end

 function run()    
    input_raw = Load_RO_weather_data(2015;stand_type=:Fertilized)
    GPP_raw = get_GPP_data(2015;stand_type=:Fertilized)
    Ec_raw = calc_Ec_data(input_raw)    
    Xₜ_raw =  Xₜ_fun(input_raw;Smin=-4.0,ΔS=17.82,τ=14.6)
    ζ = 1.19

    input_growth = input_raw[125:271]
    GPP_growth = GPP_raw[125:271]
    Ec_growth = Ec_raw[125:271]
    Xₜ_growth = Xₜ_raw[125:271]    
    
   
    GPP_model = Real[]
    Ec_model = Real[]
    GPP_old_model = Real[]
    Ec_old_model = Real[]

    for i = 1:7:147
        Xₜ =Xₜ_growth[i:i+6]
        input_week = input_growth[i:i+6]        

        optval_weekly,output_weekly,optval_day_vec,output_day_vec = run_week(input_week,Xₜ)
        
        push!(GPP_model,output_day_vec[1].GPP,
        output_day_vec[1].GPP,
        output_day_vec[2].GPP,
        output_day_vec[3].GPP,
        output_day_vec[4].GPP,
        output_day_vec[5].GPP,
        output_day_vec[6].GPP)

        push!(Ec_model,output_day_vec[1].Ec,
        output_day_vec[2].Ec,
        output_day_vec[3].Ec,
        output_day_vec[4].Ec,
        output_day_vec[5].Ec,
        output_day_vec[6].Ec,
        output_day_vec[7].Ec)

        GPP,Ec = run_week_old(input_week,Xₜ)

        push!(GPP_old_model,GPP,
        GPP,
        GPP,
        GPP,
        GPP,
        GPP,
        GPP)

        push!(Ec_old_model,Ec,
        Ec,
        Ec,
        Ec,
        Ec,
        Ec,
        Ec)
    end

    plot([day.date for day in input_growth],GPP_growth,linecolor=:blue,label="Data",xlabel="Date",ylabel="GPP")
    plot!([day.date for day in input_growth],GPP_model*ζ,linecolor=:Green,label="New Model")
    pl1 = plot!([day.date for day in input_growth],GPP_old_model*ζ,linecolor=:red,label="Old Model")

    plot([day.date for day in input_growth],Ec_growth,linecolor=:blue,xlabel="Date",ylabel="Ec")
    plot!([day.date for day in input_growth],Ec_model,linecolor=:Green)
    pl2 = plot!([day.date for day in input_growth],Ec_old_model,linecolor=:red,legends=false)

    plot(pl1,pl2,layout=(2,1))
 end

 run()
