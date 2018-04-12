# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBeam.jl/blob/master/LICENSE

""" Beam implementation for JuliaFEM. """
module FEMBeam

using FEMBase

import FEMBase: get_unknown_field_name,
                get_formulation_type,
                assemble_elements!

type Beam <: FieldProblem
end

function get_unknown_field_name(::Problem{Beam})
    return "displacement"
end

function assemble_elements!{B}(::Problem{Beam}, ::Assembly,
                               elements::Vector{Element{B}}, ::Float64)

    for element in elements
        info("Not doing anything useful right now.")
    end

end

export Beam

end
