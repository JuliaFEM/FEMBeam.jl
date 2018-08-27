# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBeam.jl/blob/master/LICENSE

"""
    FEMBeam - Beam implementation for JuliaFEM

# Supported beams
- Euler-Bernoulli beam in 3d

"""
module FEMBeam

using FEMBase, LinearAlgebra, SparseArrays

include("get_beam_stiffness_matrix_2d.jl")
include("get_beam_forces_vector_2d.jl")
include("get_beam_mass_matrix_2d.jl")
include("beam2d.jl")
export Beam2D

include("beam3d.jl")
export Beam

end
