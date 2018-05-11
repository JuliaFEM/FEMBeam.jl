# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBeam.jl/blob/master/LICENSE

using FEMBeam
using FEMBase.Test

# ABAQUS code to test:
# *NODE
#  1, 0.0, 0.0, 0.0
#  2, 0.0, 1.0, 0.0
# *ELEMENT, TYPE=B33, ELSET=EALL
#  1, 1,2
# *BEAM GENERAL SECTION, ELSET=EALL, DENSITY=100.0, SECTION=GENERAL
#  2.5E-2, 1.0E-4, 0.0E-4, 1.5E-4, 3.0E-5
#  1.0, 0.0, 0.0
#  100.0E2, 80.0E2
# *STEP
# *MATRIX GENERATE, STIFFNESS, MASS, LOAD
# *MATRIX OUTPUT, STIFFNESS, MASS, LOAD
# *DLOAD
#  1, P1, 0.0E3
#  1, P2, 0.0E3
# ** 1, PX, 54.0E3
#  1, PY, 72.0E3
# ** 1, PZ, 108.0E3
# *END STEP
# *STEP
# *STATIC
# *BOUNDARY
#  1, 1,6
# *DLOAD
#  1, P1, 0.0E3
#  1, P2, 0.0E3
# ** 1, PX, 54.0E3
#  1, PY, 72.0E3
# ** 1, PZ, 108.0E3
# *NODE PRINT
#  U, COORD
# *END STEP

L = 1.0
X = Dict(
    1 => [0.0, 0.0, 0.0],
    2 => [L, 0.0, 0.0])
beam = Element(Seg2, [1, 2])
E = 210.0e6
A = 20.0e-2
G = 84.0e6
Ix = J = 5.0e-5
I2 = Iy = 10.0e-5
I1 = Iz = 20.0e-5
update!(beam, "geometry", X)
update!(beam, "youngs modulus", E)
update!(beam, "shear modulus", G)
update!(beam, "cross-section area", A)
update!(beam, "torsional moment of inertia 1", I1)
update!(beam, "torsional moment of inertia 2", I2)
update!(beam, "polar moment of inertia", J)
update!(beam, "orientation", eye(3))
problem = Problem(Beam, "beam", 6)
add_elements!(problem, [beam])
assemble!(problem, 0.0)
K = full(problem.assembly.K)
K_expected = [
    E*A/L 0 0 0 0 0 -E*A/L 0 0 0 0 0
    0 12*E*Iz/L^3 0 0 0 6*E*Iz/L^2 0 -12*E*Iz/L^3 0 0 0 6*E*Iz/L^2
    0 0 12*E*Iy/L^3 0 -6*E*Iy/L^2 0 0 0 -12*E*Iy/L^3 0 -6*E*Iy/L^2 0
    0 0 0 G*Ix/L 0 0 0 0 0 -G*Ix/L 0 0
    0 0 -6*E*Iy/L^2 0 4*E*Iy/L 0 0 0 6*E*Iy/L^2 0 2*E*Iy/L 0
    0 6*E*Iz/L^2 0 0 0 4*E*Iz/L 0 -6*E*Iz/L^2 0 0 0 2*E*Iz/L
    -E*A/L 0 0 0 0 0 E*A/L 0 0 0 0 0
    0 -12*E*Iz/L^3 0 0 0 -6*E*Iz/L^2 0 12*E*Iz/L^3 0 0 0 -6*E*Iz/L^2
    0 0 -12*E*Iy/L^3 0 6*E*Iy/L^2 0 0 0 12*E*Iy/L^3 0 6*E*Iy/L^2 0
    0 0 0 -G*Ix/L 0 0 0 0 0 G*J/L 0 0
    0 0 -6*E*Iy/L^2 0 2*E*Iy/L 0 0 0 6*E*Iy/L^2 0 4*E*Iy/L 0
    0 6*E*Iz/L^2 0 0 0 2*E*Iz/L 0 -6*E*Iz/L^2 0 0 0 4*E*Iz/L
    ]

if !isapprox(K, K_expected)
    K_diff = K-K_expected
    K_diff[abs.(K_diff) .< 1.0e-9] = 0
    println("K")
    display(K)
    println("K_expected")
    display(K_expected)
    println("difference")
    display(K_diff)
end
@test isapprox(K, K_expected)

function solution(X2, n, q1, q2, qx, qy, qz; T=nothing)
    X = Dict(1 => [0.0, 0.0, 0.0], 2 => X2)
    beam = Element(Seg2, [1, 2])
    update!(beam, "geometry", X)
    update!(beam, "normal", n)
    update!(beam, "youngs modulus", 100.0e2)
    update!(beam, "shear modulus", 80.0e2)
    update!(beam, "cross-section area", 3.0e-2)
    update!(beam, "torsional moment of inertia 1", 1.0e-4)
    update!(beam, "torsional moment of inertia 2", 1.5e-4)
    update!(beam, "polar moment of inertia", 3.0e-5)
    update!(beam, "distributed load 1", q1)
    update!(beam, "distributed load 2", q2)
    update!(beam, "distributed load x", qx)
    update!(beam, "distributed load y", qy)
    update!(beam, "distributed load z", qz)
    if T != nothing
        update!(beam, "orientation", T)
    end
    time = 0.0
    X1, X2 = beam("geometry", time)
    L = norm(X2-X1)
    n1 = beam("normal", time)
    t = (X2-X1)/L
    n2 = cross(t,n1)
    #println("t = $t, n1  = $n1, n2 = $n2")
    problem = Problem(Beam, "example 1", 6)
    add_elements!(problem, [beam])
    assemble!(problem, time)
    K = full(problem.assembly.K)
    f = full(problem.assembly.f)
    u = zeros(12)
    fd = 7:12
    u[fd] = K[fd,fd]\f[fd]
    return K, u, f
end

include("read_abaqus.jl")

K1, u1, f1 = solution([1.0, 0.0, 0.0], [0.0, 0.0, -1.0], 0.0, 0.0, 0.0, 72.0e3, 0.0)
K1_expected = read_mtx_from_str("""
-1,1, -1,1,  3.333333333333317e+01
-1,2, -1,1, -8.333333333333492e+00
1,1, -1,1,  2.500000000000002e+01
2,1, -1,1, -2.500000000000002e+01
-1,2, -1,2,  3.333333333333316e+01
1,1, -1,2,  2.500000000000003e+01
2,1, -1,2, -2.500000000000003e+01
1,1, 1,1,  3.000000000000020e+02
2,1, 1,1, -3.000000000000020e+02
1,2, 1,2,  1.199999999999999e+01
1,6, 1,2,  5.999999999999996e+00
2,2, 1,2, -1.199999999999999e+01
2,6, 1,2,  5.999999999999996e+00
1,3, 1,3,  1.799999999999999e+01
1,5, 1,3, -8.999999999999996e+00
2,3, 1,3, -1.799999999999999e+01
2,5, 1,3, -8.999999999999996e+00
1,4, 1,4,  2.400000000000012e-01
2,4, 1,4, -2.400000000000012e-01
1,5, 1,5,  6.000000000000005e+00
2,3, 1,5,  8.999999999999996e+00
2,5, 1,5,  2.999999999999990e+00
1,6, 1,6,  4.000000000000004e+00
2,2, 1,6, -5.999999999999996e+00
2,6, 1,6,  1.999999999999993e+00
2,1, 2,1,  3.000000000000020e+02
2,2, 2,2,  1.199999999999999e+01
2,6, 2,2, -5.999999999999996e+00
2,3, 2,3,  1.799999999999999e+01
2,5, 2,3,  8.999999999999996e+00
2,4, 2,4,  2.400000000000012e-01
2,5, 2,5,  6.000000000000006e+00
2,6, 2,6,  4.000000000000004e+00""")
u1_expected = [0.000, 9000.0, 0.000, 0.000, 0.000, 1.2000E+04]
isapprox(K1, K1_expected) || begin
    println("K1")
    display(K1)
    println("K1_expected")
    display(K1_expected)
    K_diff = K1 - K1_expected
    K_diff[abs.(K_diff) .< 1.0e-9] = 0.0
    println("Difference")
    display(K_diff)
end
@test isapprox(K1, K1_expected; rtol=1.0e-9)

isapprox(u1[7:12], u1_expected) || begin
    println("f1 = $(f1[7:12])")
    println("u1 = $(u1[7:12])")
    println("u1_expected = $u1_expected")
end
@test isapprox(u1[7:12], u1_expected)

K2, u2, f2 = solution([1.0, 0.0, 0.0], [0.0, 0.0, -1.0], 0.0, 0.0, 0.0, 0.0, 108.0e3)
u2_expected = [0.000, 0.000, 9000.0, 0.000, -1.2000E+04, 0.000]
isapprox(u2[7:12], u2_expected) || begin
    println("f2 = $(f2[7:12])")
    println("u2 = $(u2[7:12])")
    println("u2_expected = $u2_expected")
end
@test isapprox(u2[7:12], u2_expected)

# without local distributed loads

K1, u1, f1 = solution([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0, 0.0, 0.0, 72.0e3, 0.0)
K1_expected = read_mtx_from_str("""
-1,1, -1,1,  3.333333333333317e+01
-1,2, -1,1, -8.333333333333492e+00
1,1, -1,1,  2.500000000000002e+01
2,1, -1,1, -2.500000000000002e+01
-1,2, -1,2,  3.333333333333316e+01
1,1, -1,2,  2.500000000000003e+01
2,1, -1,2, -2.500000000000003e+01
1,1, 1,1,  3.000000000000020e+02
2,1, 1,1, -3.000000000000020e+02
1,2, 1,2,  1.799999999999999e+01
1,6, 1,2,  8.999999999999996e+00
2,2, 1,2, -1.799999999999999e+01
2,6, 1,2,  8.999999999999996e+00
1,3, 1,3,  1.199999999999999e+01
1,5, 1,3, -5.999999999999996e+00
2,3, 1,3, -1.199999999999999e+01
2,5, 1,3, -5.999999999999996e+00
1,4, 1,4,  2.400000000000012e-01
2,4, 1,4, -2.400000000000012e-01
1,5, 1,5,  4.000000000000004e+00
2,3, 1,5,  5.999999999999996e+00
2,5, 1,5,  1.999999999999993e+00
1,6, 1,6,  6.000000000000005e+00
2,2, 1,6, -8.999999999999996e+00
2,6, 1,6,  2.999999999999990e+00
2,1, 2,1,  3.000000000000020e+02
2,2, 2,2,  1.799999999999999e+01
2,6, 2,2, -8.999999999999996e+00
2,3, 2,3,  1.199999999999999e+01
2,5, 2,3,  5.999999999999996e+00
2,4, 2,4,  2.400000000000012e-01
2,5, 2,5,  4.000000000000004e+00
2,6, 2,6,  6.000000000000006e+00""")
isapprox(K1, K1_expected) || begin
    println("K1")
    display(K1)
    println("K1_expected")
    display(K1_expected)
    K_diff = K1 - K1_expected
    K_diff[abs.(K_diff) .< 1.0e-9] = 0.0
    println("Difference")
    display(K_diff)
end
@test isapprox(K1, K1_expected; rtol=1.0e-9)

#u1_expected = [108.0, 6000.0, 1.3500E+04, 0.000, -1.8000E+04, 8000.0]
u1_expected = [0.0, 6000.0, 0.0, 0.0, 0.0, 8000.0]
isapprox(u1[7:12], u1_expected) || begin
    println("f1 = $(f1[7:12])")
    println("u1 = $(u1[7:12])")
    println("u1_expected = $u1_expected")
end
@test isapprox(u1[7:12], u1_expected)

K2, u2, f2 = solution([1.0, 0.0, 0.0], [0.0, 0.0, -1.0], 0.0, 0.0, 0*54.0e3, 72.0e3, 0*108.0e3)
K3, u3, f3 = solution([0.0, 1.0, 0.0], [1.0, 0.0, 0.0], 0.0, 0.0, 0*54.0e3, 0*72.0e3, 108.0e3)
# u4 = solution([0.0, 1.0, 0.0], [0.0, 0.0, 1.0], 0.0, 0.0, 54.0e3, 72.0e3, 108.0e3)
# u5 = solution([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 0.0, 0.0, 54.0e3, 72.0e3, 108.0e3)
# u6 = solution([0.0, 0.0, 1.0], [0.0, 1.0, 0.0], 0.0, 0.0, 54.0e3, 72.0e3, 108.0e3)

# local distributed loads
# u7 = solution([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 18.0e3, 36.0e3, 0.0, 0.0, 0.0)
# u8 = solution([1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 18.0e3, 36.0e3, 0.0, 0.0, 0.0)
# u9 = solution([0.0, 1.0, 0.0], [1.0, 0.0, 0.0], 18.0e3, 36.0e3, 0.0, 0.0, 0.0)
# u10 = solution([0.0, 1.0, 0.0], [0.0, 0.0, 1.0], 18.0e3, 36.0e3, 0.0, 0.0, 0.0)
# u11 = solution([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 18.0e3, 36.0e3, 0.0, 0.0, 0.0)
# u12 = solution([0.0, 0.0, 1.0], [0.0, 1.0, 0.0], 18.0e3, 36.0e3, 0.0, 0.0, 0.0)


#u2_expected = [108.0, 9000.0, 9000.0, 0.000, -1.2000E+04, 1.2000E+04]
#u3_expected = [4500.0, 144.0, 1.3500E+04, 1.8000E+04, 0.000, -6000.0]

#@test isapprox(u2, u2_expected)
#@test isapprox(u3, u3_expected)
