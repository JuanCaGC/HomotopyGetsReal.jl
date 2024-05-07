module edgeModule

include("cell.jl")
using .cellModule

struct Edge <: Cell
    dimension::Int
    function Edge(dimension::Int=1)
        new(dimension)
    end

end
end