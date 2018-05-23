# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBeam.jl/blob/master/LICENSE

""" Beam implementation for JuliaFEM. """
module FEMBeam

using Reexport
@reexport using FEMBase

include("get_beam_stiffness_matrix_2d.jl")
include("get_beam_forces_vector_2d.jl")
include("get_beam_mass_matrix_2d.jl")

import FEMBase: get_unknown_field_name,
                get_formulation_type,
                assemble_elements!

type Beam2D <: FieldProblem
end

function get_unknown_field_name(::Problem{Beam2D})
    return "displacement"
end

function assemble_elements!{B}(::Problem{Beam2D}, ::Assembly,
                               elements::Vector{Element{B}}, ::Float64)

    for element in elements
        info("Not doing anything useful right now.")
    end

end

export Beam2D

include("beam3d.jl")

export Beam

end
