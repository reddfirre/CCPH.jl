#Bisection method for finding one root to f in the interval [a,b]
function bisection(f::Function, a::AbstractFloat, b::AbstractFloat;
    tol::AbstractFloat=1e-5, maxiter::Integer=200)
    fa = f(a)
    fa*f(b) <= 0 || error("No real root in [a,b]")
    i = 0
    local c
    while b-a > tol
    i += 1
    i != maxiter || error("Max iteration exceeded")
    c = (a+b)/2
    fc = f(c)
    if fc == 0
    break
    elseif fa*fc > 0
    a = c  # Root is in the right half of [a,b].
    fa = fc
    else
    b = c  # Root is in the left half of [a,b].
    end
    end
    return c
end

#Fixed-point iteration method for finding a solution to the eqution f(x) = x
function fixpoint(f::Function,x0::Real;Max_step::Integer=200,Abs_tol::Real=10^-4)
    x_old = x0
    x_new = f(x_old)
    abs_error = abs(x_old-x_new)  
    step_nr = 1 
    while abs_error>Abs_tol
        step_nr < Max_step||error("fixpoint: |f(x)-x|>Abs_tol after Max_step iterations")
        x_old = x_new       
        x_new = f(x_old)
        abs_error = abs(x_old-x_new)
        step_nr += 1        
    end
    return x_new
end

#Function for finding the real roots for the real polynomial f(x) = ax³+bx²+cx+d
function solvecubic(a::Real, b::Real, c::Real, d::Real)
    function one_real_root(a::Real, b::Real, g::Real, h::Real)
        R = -g/2 + sqrt(h)
        S = cbrt(R)
        T = -g/2 - sqrt(h)
        U = cbrt(T)
        x1 = (S + U) - (b/3a)
        x2 = -(S + U)/2 - (b/3a) - im*(S-U)*sqrt(3) / 2
        x3 = -(S + U)/2 - (b/3a) + im*(S-U)*sqrt(3) / 2
        return Union{Real,Complex}[x1, x2, x3]
    end  
    
    function three_real_root(a::Real, b::Real, g::Real, h::Real)
        i = (g^2 / 4 - h)^(1/2)
        j = cbrt(i)
        K = acos(-(g/2i))
        L = -j
        M = cos(K/3)
        N = sqrt(3) * sin(K/3)
        P = -(b/3a)
        x1 = 2j * cos(K/3) - (b/3a)
        x2 = L * (M + N) + P
        x3 = L * (M - N) + P
        return sort([x1, x2, x3])
    end

    if a == zero(typeof(a))
        return Base.error("Can't solve when a=0 in a*x³ + b*x² + c*x + d =0")
    end
    f = (3c/a - (b/a)^2) / 3
    g = (2(b/a)^3 - (9*b*c/a^2) + 27d/a) / 27
    h = g^2 / 4 + f^3 / 27
    if f == zero(typeof(f)) && g == zero(typeof(g)) && h == zero(typeof(h))
        x1 = x2 = x3 = -cbrt(d/a)
        return [x1, x2, x3]
    elseif  h <= zero(typeof(h))
        return three_real_root(a, b, g, h)
    elseif h > zero(typeof(h))
        return one_real_root(a, b, g, h)
    end

    return Base.error("Can't solve when a=$(a) b=$(b) c=$(c) d=$(d)")
end

#Dynamic viscosity of water (Pa s) at T degree celsius
#Equation 5 in Korson et al. (1969)
ηₜ(T::AbstractFloat) = 1.0020*10^((1.1709*(20-T)-0.001827*(T-20)^2)/(T+89.93)-3)

#Saturation Vapor Pressure (Pa) at T degree celsius
#Equation 21 in Alduchov and Eskridge (1996)
SVPₜ(T::AbstractFloat) = 610.94*exp(17.625*T/(243.04+T))

#Vapour-pressure deficit (Pa) at T degree celsius and RH relative humidity (%)
VPDₜ(T::AbstractFloat,RH::AbstractFloat) = SVPₜ(T)*(1-RH/100)

#Calcualte day-light hours
function daylighthour(L::Real,d::Real)
    #Jenkins, 2013, The Sun's position in the sky (eq 17)
    #L = latitude (rad)
    #d = day of the year (d=1 is january 1)
    M = -0.0410+0.017202*d
    Φ = -1.3411+M+0.0334*sin(M)+0.0003*sin(2*M)
    ϵ = sin(23.4*pi/180)
    Φs = sin(Φ)
    return 24*(1-1/pi*acos(tan(L)*ϵ*Φs/sqrt(1-ϵ^2*Φs^2)))
end