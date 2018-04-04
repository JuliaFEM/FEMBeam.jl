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
    ne=2
    BCs=[1 , 2 , 3]
    qt=0
    qn=-5.0e3
    F=-0.0e3

    U_expected = [0.0, 1.15491e-17, -3.46474e-17, 0.0, -0.0337302, -0.0555556, 0.0, -0.0952381, -0.0634921]

    K_expected =
    [2.94e8 0.0      0.0     -2.94e8  0.0         0.0         0.0    0.0      0.0;
     0.0    1.26e6   630000.0 0.0    -1.26e6      630000.0    0.0    0.0      0.0;
     0.0    630000.0 420000.0 0.0    -630000.0    210000.0    0.0    0.0      0.0;
    -2.94e8 0.0      0.0      5.88e8  0.0         0.0        -2.94e8 0.0      0.0;
     0.0   -1.26e6  -630000.0 0.0     2.52e6      2.32831e-10 0.0   -1.26e6   630000.0;
     0.0    630000.0 210000.0 0.0     2.32831e-10 840000.0    0.0   -630000.0 210000.0;
     0.0    0.0      0.0     -2.94e8  0.0         0.0         2.94e8 0.0      0.0;
     0.0    0.0      0.0      0.0    -1.26e6     -630000.0    0.0    1.26e6  -630000.0;
     0.0    0.0      0.0      0.0     630000.0    210000.0    0.0   -630000.0 420000.0]

     la_expected=[0.0, -10000.0, -10000.0]

    U,K,la = FEMBeam.fembeam(l,E,I,A,ro,ne,BCs,qt,qn,F)

    @test isapprox(U, U_expected)
    @test isapprox(K, K_expected)
    @test isapprox(la, la_expected)
end
