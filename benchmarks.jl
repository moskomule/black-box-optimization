# BBOB benchmarks
# http://coco.lri.fr/downloads/download15.03/bbobdocfunctions.pdf

module Benchmark
using LinearAlgebra

export rotation_matrix
function rotation_matrix(D)
    return qr(rand(D, D)).Q
end
    
function _t_asy(x, beta)
    new_x = copy(x)
    D = length(x)
    for i in 1:D
        if x[i] > 0
            new_x[i] = x[i]^(1 + beta * (i - 1) / (D - 1) * sqrt(x[i]))
        end
    end
    return new_x
end

function _t_osz(x)
    new_x = zero(x)
    D = length(x)
    for i in 1:D
        c1 = x[i] > 0 ? 10 : 5.5
        c2 = x[i] > 0 ? 7.9 : 3.1
        _x = x[i] != 0 ? log(abs(x[i])) :  0
        new_x[i] = sign(x[i]) * exp(_x + 0.049 * (sin(c1 * _x) + sin(c2 * _x)))
    end
    return new_x
end

function _lambda(D, alpha)
    ret = zeros(D)
    for i in 1:D
        ret[i] = alpha^(0.5 * (i - 1) / (D - 1))
    end
    return ret
end

function _f_pen(x)
    new_x = zero(x)
    for i in 1:length(x)
        t = abs(x[i]) - 5
        if t > 0
            new_x[i] = t^2
        end
    end
    return sum(new_x)

end

export sphere
function sphere(x, x0)
    return sum((x - x0).^2)
end

export ellipsoidal
function ellipsoidal(x, x0)
   # x_opt = 0
    f = 0
    D = length(x)
    z = _t_osz(x - x0)
    for i in 1:D
        f += 10^((i - 1) / (D - 1)) * z[i]^2
    end
    return f
end

export rastrigin
function rastrigin(x, x0)
    # x_opt = 0
    D = length(x)
    z = _lambda(D, 10.) .* _t_asy(_t_osz(x - x0), 0.2)
    return 10 * (D - sum(cos.(2pi * z))) + sum(z.^2)
end

export buche_rastigin
function buche_rastigin(x, x0)
    # x_opt = 0
    s = zero(x)
    D = length(x)
    for i in 1:D
        if x[i] > 0 & i % 2 == 1
            s[i] = 10 * 10^(0.5 * (i - 1) / (D - 1))
        else
            s[i] = 10^(0.5 * (i - 1) / (D - 1))
        end
    end

    z = s .* _t_osz(x - x0)

    return 10 * (D - sum(cos.(2pi * z))) + sum(z.^2) + 100 * _f_pen(x)
    
end

export linear_slope
function linear_slope(x, x0)
    s = zero(x) 
    D = length(x)
    for i in 1:D
        s[i] = sign(x0[i] * 10^((i - 1) / (D - 1)))
    end
    z = copy(x0)
    for i in 1:D
        if x0[i] * x[i] < 25
            z[i]  = x[i]
        end
    end
    return sum(5 * abs.(s) - s .* z)
end

export attractive_sector
function attractive_sector(x, x0, Q)
    D = length(x)
    z = Q * diagm(_lambda(D, 10.)) * Q * (x - x0)
    s = zero(x) .+ 1
    for i in 1:D
        if z[i] * x0[i] > 0
            s[i] = 100
        end
    end
    return _t_osz([sum((s .* z).^2)])[1]^0.9
end

export step_elipsoidal
function step_elipsoidal(x, x0, R1, R2)
    D = length(x)
    z_hat = diagm(_lambda(D, 10.)) * R1 * (x - x0)
    z_til = zero(x)
    for i in 1:D
        _hat = z_hat[i]
        z_til[i] = _hat > 0.5 ? floor(0.5 + _hat) : floor(0.5 + _hat) / 10
    end
    z = R2 * z_til
    f = 0
    for i in 1:D
        f += 10^(2 * (i - 1) / (D - 1)) * z[i]^2
    end
    return 0.1 * max(z_hat[1] / 10^4, f) + _f_pen(x)
end

export rosenbrock
function rosenbrock(x, x0)
    D = length(x)
    z = max(1, sqrt(D) / 8) * (x - x0) .+ 1.0
    f = 0
    for i in 1:D - 1
        f += (100 * (z[i] - z[i + 1])^2 + (z[i] - 1)^2) 
    end
    return f
end


export rosenbroch_rotated
function rosenbroch_rotated(x, x0, R)
    D = length(x)
    z = max(1, sqrt(D) / 8) * R * (x - x0) .+ 0.5
    f = 0
    for i in 1:D - 1
        f += (100 * (z[i] - z[i + 1])^2 + (z[i] - 1)^2) 
    end
    return f
end

export ellipsoidal_rotated
function ellipsoidal_rotated(x, x0, R)
    z = _t_osz(R * (x - x0))
    D = length(x)
    f = 0
    for i in 1:D
        f += 10^(6 * (i - 1) / (D - 1)) * z[i]^2
    end
    return f
end

export discus
function discus(x, x0, R)
    z = _t_osz(R * (x - x0))
    return 10^6 * z[1]^2 + sum(z[2:length(x)].^2)
end

export bent_cigar
function bent_cigar(x, x0, R)
    z = R * _t_asy(R(x - x0), 0.5)
    return z[1]^2 + 10^6 * sum(z[2:length(x)].^2)
end


export sharp_ridge
function sharp_ridge(x, x0, R, Q)
    D = length(x)
    z = R * diagm(_lambda(D, 10.)) * Q * (x - x0)
    
    return z[1]^2 + 100 * sqrt(sum(z[2:D].^2))
end

export different_powers
function different_powers(x, x0, R)
    z = R * (x - x0)
    D = length(x)
    f = 0
    for i in 1:D
        f += abs(z[i])^(2 + 4 * (i - 1) / (D - 1))
    end
    return sqrt(f)
end

export rastrigin_rotated
function rastrigin_rotated(x, x0, R, Q)
    D = length(x)
    z = R * diagm(_lambda(D, 10.)) * Q * _t_asy(_t_osz(x - x0), 0.2)
    return 10 * (D - sum(cos.(2pi * z))) + sum(z.^2)
end

export weierstrass
function weierstrass(x, x0, R, Q)
    D = length(x)
    z = R * diagm(_lambda(D, 1 / 100)) * Q * _t_osz(R * (x - x0))
    f_0 = 0
    for k in 0:11
        f_0 += 2^(-k) * cos(2pi * 3^k / 3)
    end
    f = 0
    for i in 1:D
        for k in 0:11
            f += 2^(-k) * cos(2pi * 3^k * (z[i] + 0.5))
        end
    end
    return 10 * (1 / D * f - f_0)^3 + 10 / D * _f_pen(x)
end

export schaffers_f7
function schaffers_f7(x, x0, R, Q)
    D = length(x)
    z = diagm(D, 10.) * Q * _t_asy(R(x - x0), 0.5)
    s = zero(x)
    for i in 1:D - 1
        s[i] = sqrt(z[i]^2 + z[i + 1]^2)
    end
    # I don't know this is correct
    s[D] = sqrt(z[0]^2 + z[D]^2)
    f = 0
    for i in 1:D - 1
        f += sqrt(s[i]) * (1 + (sin(50 * s[i]^0.2)^2))
    end
    return (f / (D - 1))^2 + 10 * _f_pen(x)
end

function griewank_rosenbrock()
    
end

function schwefel_226()
    
end

function gallaghers_gaussian_101me_peaks()
    
end

function gallaghers_gaussian_21hi_peaks()
    
end

function katsuura()
    
end

function lunacek_bi_rastrigin()
    
end

# module
end