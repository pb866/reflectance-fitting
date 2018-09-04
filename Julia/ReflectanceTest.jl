using ALS
using Test

# Unit test for Reflectance
@testset "Reflectance" begin
    th1 = [t for t in 1.0:5.0]
    th2 = th1 .+ 0.1
    rf1 = cos.(th1)
    rf2 = cos.(th2)
    refl1 = Reflectance(th1, rf1, 25.0);
    refl2 = Reflectance(th2, rf2, 25.0);
    refls = refl1+refl2
    @test length(refl1.theta)+length(refl2.theta) == length(refls.theta)
    ok = true
    lastth = -1000.0
    for th in refls.theta
        if th < lastth
            ok = false
        end
    end
    @test ok
end
