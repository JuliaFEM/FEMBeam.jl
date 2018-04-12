# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBeam.jl/blob/master/LICENSE

using FEMBase
using FEMBeam

using Base.Test

@testset "FEMBeam test" begin
X1=[0.0,0.0]; X2=[6.0,0.0]
E=210.0e9
I=50/100^4
A=14/100^2
ro=7800
qt=-5.0e3
qn=0
f=zeros(6,1)

k_expected =
[4.9e7 0.0 0.0 -4.9e7 0.0 0.0;
0.0 5833.33 17500.0 0.0 -5833.33 17500.0;
0.0 17500.0 70000.0 0.0 -17500.0 35000.0;
-4.9e7 0.0 0.0 4.9e7 0.0 0.0;
0.0 -5833.33 -17500.0 0.0 5833.33 -17500.0;
0.0 17500.0 35000.0 0.0 -17500.0 70000.0]

 m_expected=
 [21.84 0.0 0.0 10.92 0.0 0.0;
 0.0 24.336 20.592 0.0 8.424 -12.168;
 0.0 20.592 22.464 0.0 12.168 -16.848;
 10.92 0.0 0.0 21.84 0.0 0.0;
 0.0 8.424 12.168 0.0 24.336 -20.592;
 0.0 -12.168 -16.848 0.0 -20.592 22.464]

 f_expected=[0.0, -15000.0, -15000.0, 0.0, -15000.0, 15000.0]

k,m,f = FEMBeam.fembeam(X1,X2,E,I,A,ro,qt,qn,f)

@test isapprox(k, k_expected)
@test isapprox(m, m_expected)
@test isapprox(f, f_expected)



k,m,f = FEMBeam.fembeam(X1,X2,E,I,A,ro,qt,qn,f)


end
