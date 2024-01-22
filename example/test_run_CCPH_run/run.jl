using CCPH

H = 19.07 #m
N = 850.0/10000 # #tree/m²
LAI = 2.45 #-

Tmean = 14.1 #°C
Tmin = 7.9 #°C
Tmax = 19.7 #°C
VP = 1320.0 #Pa
Radₜₒₜ = 10.3*10^6 #J m⁻²
I₀ₜₒₜ = Radₜₒₜ*2.3*10^-6 #mol m⁻²
θₛ = 0.256 #-
daylength = CCPH.daylighthour(64*pi/180,217)*3600 #Seconds

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
treepar = TreePar(Nₛ=0.013,α_max=0.14,a_Jmax=0.011,b_Jmax=0.0)
treesize = TreeSize(H,LAI,N)
hydPar = HydraulicsPar(Kₓₗ₀=0.00054)
kinetic = PhotoKineticRates()
photo = PhotoPar(kinetic,env.Tₐ)  

model = CCPHStruct(cons,env,treepar,treesize,photo,hydPar)

gₛ₁ = 0.13 
gₛ₂ = 0.13
Nₘ_f = 0.02

@show modeloutput = CCPH_run!(gₛ₁,gₛ₂,Nₘ_f,daylength,kinetic,envfun,model)

@show gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt = CCPHOpt(daylength,kinetic,envfun,model)

@show CCPH_run!(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt,daylength,kinetic,envfun,model)

@show t₁ = daylength*asin(2/π)/(2*π)
@show t₂ = daylength/4+t₁

@show VP/CCPH.SVPₜ(Tₐ_fun(t₁))
@show VP/CCPH.SVPₜ(Tₐ_fun(t₂))
@show VPD_fun(t₁)
@show VPD_fun(t₂)