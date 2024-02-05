using CCPH

function old_obj_fun(x,model::CCPH.CCPHStruct)
    gₛ,Nₘ_f = x
    output = CCPH.CCPH_inst_run(gₛ,Nₘ_f,model)
    return -output.Gain
end

H = 19.07 #m
N = 850.0/10000 # #tree/m²
LAI = 2.45 #-

Tmean = 11.2 #°C
Tmin = 4.7 #°C
Tmax = 17.7 #°C
VP_data = 1110 #Pa
@show VPmin = CCPH.SVPₜ(Tmin) #Pa
@show VP = min(VP_data,VPmin)
Radₜₒₜ = 10.4*10^6 #J m⁻²
I₀ₜₒₜ = Radₜₒₜ*2.3*10^-6 #mol m⁻²
θₛ = 0.2569903364#-
daylength = CCPH.daylighthour(64*pi/180,199)*3600 #Seconds

@show ψₛ = CCPH.θₛ2ψₛ(θₛ) #MPa
@show Cₐ = 400.0/10.0 #Pa
@show P = 1.0*10^5 #Pa
@show Tₐ = Tmean #°C
@show VPD =  CCPH.VPDₜ(Tₐ,VP) #Pa
@show I₀ = I₀ₜₒₜ/daylength #mol m⁻² s⁻¹

ψₛ_fun(t) = ψₛ
Cₐ_fun(t) = Cₐ
P_fun(t) = P
Tₐ_fun(t) = CCPH.Tₐ_fun(t,Tmin,Tmax,daylength)
VPD_fun(t) = CCPH.VPDₜ(Tₐ_fun(t),VP)
I₀_fun(t) = CCPH.I₀_fun(t,I₀ₜₒₜ,daylength)

cons = Constants()
envfun = EnvironmentFunStruct(I₀_fun,Cₐ_fun,P_fun,Tₐ_fun,VPD_fun,ψₛ_fun)
env = EnvironmentStruct(I₀,Cₐ,P,Tₐ,VPD,ψₛ)
treepar = TreePar(Xₜ=0.9,Nₛ=0.013,α_max=0.14,a_Jmax=0.011,b_Jmax=0.0)
treesize = TreeSize(H,LAI,N)
hydPar = HydraulicsPar(Kₓₗ₀=0.00054)
kinetic = PhotoKineticRates()
photo = PhotoPar(kinetic,env.Tₐ)  

model = CCPHStruct(cons,env,treepar,treesize,photo,hydPar)

#--Old Model--
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
@show Ec

GPP = output.GPP*daylength
GPP *= model.cons.M_C*1000
@show GPP

#--New Model--
gₛ₁_lim_hi,gₛ₂_lim_hi = CCPH.SDM2_get_gₛ_lim!(daylength,model,kinetic,envfun)
gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt = CCPH.CCPHOpt(daylength,kinetic,envfun,model;gₛ₁_lim_hi=gₛ₁_lim_hi,gₛ₂_lim_hi=gₛ₂_lim_hi)
output_new = CCPH_run!(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt,daylength,kinetic,envfun,model)

#Output from new model
@show output_new.Ec
@show output_new.GPP