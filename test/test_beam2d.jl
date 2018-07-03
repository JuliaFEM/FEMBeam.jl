# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBeam.jl/blob/master/LICENSE

using FEMBeam
using FEMBase.Test
beam_element = Element(Seg2, [1, 2])
problem = Problem(Beam2D, "test problem", 3)
add_elements!(problem, [beam_element])

@test get_unknown_field_name(problem) == "displacement"
@test assemble!(problem, 0.0) == nothing
