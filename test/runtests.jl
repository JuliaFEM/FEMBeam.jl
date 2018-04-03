# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBeam.jl/blob/master/LICENSE

using FEMBase
using FEMBeam

using Base.Test

@testset "FEMBeam test" begin
l=2.0
E=210.0e9
I=50/100^4
A=14/100^2
ro=1
ne=1
BCs=[1 , 2 , 3]
qt=5e3
qn=-5.0e3
F=-0.0e3

    U_expected = [0.0, 1.44364e-17, 8.66185e-18, 3.40136e-5, -0.0952381, -0.0634921, 10000.0, -10000.0, -10000.0]
    K_expected = [1.47e8 0.0 0.0 -1.47e8 0.0 0.0 1.0 0.0 0.0;
    0.0 157500.0 157500.0 0.0 -157500.0 157500.0 0.0 1.0 0.0;
    0.0 157500.0 210000.0 0.0 -157500.0 105000.0 0.0 0.0 1.0;
    -1.47e8 0.0 0.0 1.47e8 0.0 0.0 0.0 0.0 0.0;
    0.0 -157500.0 -157500.0 0.0 157500.0 -157500.0 0.0 0.0 0.0;
    0.0 157500.0 105000.0 0.0 -157500.0 210000.0 0.0 0.0 0.0;
    1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0]
    # M_expected = [1,1,1]

    U,K = FEMBeam.fembeam(l,E,I,A,ro,ne,BCs,qt,qn,F)

    @test isapprox(K, K_expected)
    @test isapprox(U, U_expected)
end
