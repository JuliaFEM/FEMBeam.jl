X = Dict(
    1 => [0.0, 0.0, 0.0],
    2 => [3.0, 0.0, 0.0],
    3 => [0.0, 0.0, -3.0],
    4 => [0.0, -4.0, 0.0])

beam1 = Element(Seg2, [1, 2])
update!(beam1, "geometry", X)
update!(beam1, "youngs modulus", 100.0e6)
update!(beam1, "shear modulus", 80.0e6)
update!(beam1, "cross-section area", 2.5e-2)
update!(beam1, "torsional moment of inertia 1", 1.0e-4)
update!(beam1, "torsional moment of inertia 2", 1.5e-4)
update!(beam1, "polar moment of inertia", 3.0e-5)
#update!(beam1, "distributed load 1", 1.0e3)

beam2 = Element(Seg2, [1, 3])
update!(beam2, "geometry", X)
update!(beam2, "youngs modulus", 120.0e6)
update!(beam2, "shear modulus", 90.0e6)
update!(beam2, "cross-section area", 1.5e-2)
update!(beam2, "torsional moment of inertia 1", 1.2e-4)
update!(beam2, "torsional moment of inertia 2", 2.0e-4)
update!(beam2, "polar moment of inertia", 4.0e-5)
#update!(beam2, "distributed load 2", 2.0e3)

beam3 = Element(Seg2, [1, 4])
update!(beam3, "geometry", X)
update!(beam3, "youngs modulus", 150.0e6)
update!(beam3, "shear modulus", 70.0e6)
update!(beam3, "cross-section area", 1.8e-2)
update!(beam3, "torsional moment of inertia 1", 1.3e-4)
update!(beam3, "torsional moment of inertia 2", 1.4e-4)
update!(beam3, "polar moment of inertia", 5.0e-5)
#update!(beam3, "distributed load x", 1.0e3)
#update!(beam3, "distributed load y", 2.0e3)
#update!(beam3, "distributed load z", 3.0e3)

load1 = Element(Poi1, [1])
update!(load1, "point force 1",  1.0e3)
update!(load1, "point force 2",  2.0e3)
update!(load1, "point force 3",  3.0e3)
update!(load1, "point moment 1", 4.0e3)
update!(load1, "point moment 2", 5.0e3)
update!(load1, "point moment 3", 6.0e3)

# update!(beam1, "normal", [0.0, 0.0, 1.0])
# update!(beam2, "normal", [0.0, 1.0, 0.0])
# update!(beam3, "normal", [0.0, 0.0, 1.0])

update!(beam1, "normal", [0.0, 1.0, 0.0])
update!(beam2, "normal", [1.0, 0.0, 0.0])
update!(beam3, "normal", [1.0, 0.0, 0.0])

problem = Problem(Beam, "example 1", 6)
add_elements!(problem, [beam1, beam2, beam3, load1])
time = 0.0
assemble!(problem, time)
K = full(problem.assembly.K)[1:6,1:6]
f = full(problem.assembly.f)[1:6]
u = K\f
println("u = $u")
println("f = $f")
#u_expected = [1.1040E-03, -5.5843E-04, 7.9695E-03, 7.5682E-02, 0.1270, 0.1762]

# Without distributed loads:
# [0.00110405, -0.00055843, 0.00796946, 0.0756821, 0.127019, 0.17617]
# ABAQUS solution:
# 1.1040E-03 -5.5843E-04  7.9695E-03  7.5682E-02  0.1270      0.1762
