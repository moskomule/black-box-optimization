# minimum test
include("benchmarks.jl")

for dim in 3:10
    x = zeros(dim)
    x0 = zeros(dim)
    Q = Benchmarks.rotation_matrix(dim)
    R = Benchmarks.rotation_matrix(dim)
    
    if Benchmarks.sphere(x, x0) != 0
        println("sphere function fails when $dim")
    end

    if Benchmarks.ellipsoidal(x, x0) != 0
        println("ellipsoidal function fails when $dim")
    end

    if Benchmarks.rastrigin(x, x0) != 0
        println("rastrigin function fails when $dim")
    end

    if Benchmarks.buche_rastigin(x, x0) != 0
        println("buche_rastigin function fails when $dim")
    end

    if Benchmarks.linear_slope(x, x0) != 0
        println("linear_slope function fails when $dim")
    end

    if Benchmarks.attractive_sector(x, x0, Q) != 0
        println("attractive_sector fails when $dim")
    end

    if Benchmarks.step_ellipsoidal(x, x0, R, Q) != 0
        println("step_ellipsoidal fails when $dim")
    end

    if Benchmarks.rosenbrock(x, x0) != 0
        println("rosenbrock fails when $dim")
    end

    if Benchmarks.ellipsoidal_rotated(x, x0, R) != 0
        println("ellipsoidal_rotated fails when $dim")
    end

    if Benchmarks.discus(x, x0, R) != 0
        println("discus fails when $dim")
    end

    if Benchmarks.bent_cigar(x, x0, R) != 0
        println("bent_cigar fails when $dim")
    end

    if Benchmarks.sharp_ridge(x, x0, R, Q) != 0
        println("sharp_ridge fails when $dim")
    end

    if Benchmarks.different_powers(x, x0, R) != 0
        println("different_powers fails when $dim")
    end

    if Benchmarks.rastrigin_rotated(x, x0, R, Q) != 0
        println("rastrigin_rotated fails when $dim")
    end

    if Benchmarks.weierstrass(x, x0, R, Q) != 0
        println("weierstrass fails when $dim")
    end

    if Benchmarks.schaffers_f7(x, x0, R, Q) != 0
        println("schaffers_f7 fails when $dim")
    end

end