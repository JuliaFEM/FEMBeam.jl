# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBeam.jl/blob/master/LICENSE

using FEMBase
using FEMBeam

using Base.Test

@testset "Beam 1 Stiffness matrix" begin
X1=[0.0,0.0]; X2=[0.0,6.5]
E=210.0e9
I=50.8e-3*101.6e-3^3/12
A=50.8e-3*101.6e-3
k = FEMBeam.get_beam_stiffness_matrix_2d(X1,X2,E,I,A)
k_expected=
[40740.3     1.02079e-8  -1.32406e5    -40740.3    -1.02079e-8  -1.32406e5;
 1.02079e-8  1.66749e8    8.10752e-12  -1.02079e-8 -1.66749e8    8.10752e-12;
-1.32406e5   8.10752e-12  573759.0      1.32406e5  -8.10752e-12  2.8688e5;
-40740.3    -1.02079e-8   1.32406e5     40740.3     1.02079e-8   1.32406e5;
-1.02079e-8 -1.66749e8   -8.10752e-12   1.02079e-8  1.66749e8   -8.10752e-12;
-1.32406e5   8.10752e-12  2.8688e5      1.32406e5  -8.10752e-12  573759.0]
@test isapprox(k, k_expected, rtol=0.0001)
end

@testset "Beam 1 force vector" begin
X1=[0.0,0.0]; X2=[0.0,6.5]
qt=-500
qn=0
f=zeros(6,1)
f = FEMBeam.get_beam_forces_vector_2d(X1,X2,qt,qn,f)
f_expected=[1625.0, -9.95026e-14, -1760.42, 1625.0, -9.95026e-14, 1760.42]
@test isapprox(f, f_expected, rtol=0.0001)
end

@testset "Beam 1 mass matrix" begin
X1=[0.0,0.0]; X2=[0.0,6.5]
A=50.8e-3*101.6e-3
ro=7800
m = FEMBeam.get_beam_mass_matrix_2d(X1,X2,A,ro)
m_expected=
[97.1943 -6.10403e-16 -89.0948 33.6442 6.10403e-16 52.6469;
-6.10403e-16 87.2256 5.45548e-15 6.10403e-16 43.6128 -3.22369e-15;
-89.0948 5.45548e-15 105.294 -52.6469 3.22369e-15 -78.9703;
33.6442 6.10403e-16 -52.6469 97.1943 -6.10403e-16 89.0948;
6.10403e-16 43.6128 3.22369e-15 -6.10403e-16 87.2256 -5.45548e-15;
52.6469 -3.22369e-15 -78.9703 89.0948 -5.45548e-15 105.294]
@test isapprox(m, m_expected, rtol=0.0001)
end

@testset "Beam 2 Stiffness matrix" begin
X1=[0.0,6.5]; X2=[8.0,6.5]
E=250.0e9
I=75.0e-3*100.0e-3^3/12
A=75.0e-3*100.0e-3
k = FEMBeam.get_beam_stiffness_matrix_2d(X1,X2,E,I,A)
k_expected=
[2.34375e8 0.0 0.0 -2.34375e8 0.0 0.0;
0.0 36621.1 1.46484e5 0.0 -36621.1 1.46484e5;
0.0 1.46484e5 781250.0 0.0 -1.46484e5 390625.0;
-2.34375e8 0.0 0.0 2.34375e8 0.0 0.0;
0.0 -36621.1 -1.46484e5 0.0 36621.1 -1.46484e5;
0.0 1.46484e5 390625.0 0.0 -1.46484e5 781250.0]
@test isapprox(k, k_expected, rtol=0.0001)
end

@testset "Beam 2 force vector" begin
X1=[0.0,6.5]; X2=[8.0,6.5]
qt=-750
qn=0
f=zeros(6,1)
f = FEMBeam.get_beam_forces_vector_2d(X1,X2,qt,qn,f)
f_expected=[0.0, -3000.0, -4000.0, 0.0, -3000.0, 4000.0]
@test isapprox(f, f_expected, rtol=0.0001)
end

@testset "Beam 2 mass matrix" begin
X1=[0.0,6.5]; X2=[8.0,6.5]
A=75e-3*100e-3
ro=10000
m = FEMBeam.get_beam_mass_matrix_2d(X1,X2,A,ro)
m_expected=
[200.0 0.0 0.0 100.0 0.0 0.0;
0.0 222.857 251.429 0.0 77.1429 -148.571;
0.0 251.429 365.714 0.0 148.571 -274.286;
100.0 0.0 0.0 200.0 0.0 0.0;
0.0 77.1429 148.571 0.0 222.857 -251.429;
0.0 -148.571 -274.286 0.0 -251.429 365.714]
@test isapprox(m, m_expected, rtol=0.0001)
end

@testset "Beam 3 Stiffness matrix" begin
X1=[8.0,6.5]; X2=[8.0,0.0]
E=160.0e9
I=50.0e-3*50.0e-3^3/12
A=50.0e-3*50.0e-3
k = FEMBeam.get_beam_stiffness_matrix_2d(X1,X2,E,I,A)
k_expected=
[3641.33 -3.76792e-9 11834.3 -3641.33 3.76792e-9 11834.3;
-3.76792e-9 6.15385e7 7.24643e-13 3.76792e-9 -6.15385e7 7.24643e-13;
11834.3 7.24643e-13 51282.1 -11834.3 -7.24643e-13 25641.0;
-3641.33 3.76792e-9 -11834.3 3641.33 -3.76792e-9 -11834.3;
3.76792e-9 -6.15385e7 -7.24643e-13 -3.76792e-9 6.15385e7 -7.24643e-13;
11834.3 7.24643e-13 25641.0 -11834.3 -7.24643e-13 51282.1]
@test isapprox(k, k_expected, rtol=0.0001)
end

@testset "Beam 3 force vector" begin
X1=[8.0,6.5]; X2=[8.0,0.0]
qt=0
qn=0
f=zeros(6,1)
f = FEMBeam.get_beam_forces_vector_2d(X1,X2,qt,qn,f)
f_expected=zeros(6,1)
@test isapprox(f, f_expected, rtol=0.0001)
end

@testset "Beam 3 mass matrix" begin
X1=[8.0,6.5]; X2=[8.0,0.0]
A=50.0e-3*50.0e-3
ro=6000
m = FEMBeam.get_beam_mass_matrix_2d(X1,X2,A,ro)
m_expected=
[36.2143 2.27434e-16 33.1964 12.5357 -2.27434e-16 -19.6161;
2.27434e-16 32.5 2.03269e-15 -2.27434e-16 16.25 -1.20114e-15;
33.1964 2.03269e-15 39.2321 19.6161 1.20114e-15 -29.4241;
12.5357 -2.27434e-16 19.6161 36.2143 2.27434e-16 -33.1964;
-2.27434e-16 16.25 1.20114e-15 2.27434e-16 32.5 -2.03269e-15;
-19.6161 -1.20114e-15 -29.4241 -33.1964 -2.03269e-15 39.2321]
@test isapprox(m, m_expected, rtol=0.0001)
end

@testset "test beam 3D stiffness" begin
    include("test_beam3d_ex1.jl")
end

@testset "test beam 3D mass matrix" begin
    include("test_beam3d_mass_matrix.jl")
end

@testset "test supports" begin
    include("test_supports.jl")
end

@testset "test_rotation_matrix.jl" begin include("test_rotation_matrix.jl") end
