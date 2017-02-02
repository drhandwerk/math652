include "Vertex.jl"

"""
  Edge(p1, p2)

Creates an edge between two vertices 'p1' and 'p2'.
"""
type Edge(p1::Vertex, p2::Vertex) <: AbstractArray{Float64, 1}
  p1::Vertex
  p2::Vertex

  Edge(p1::Vertex, p2::Vertex) = new(p1, p2)
end
Base.size(E::Edge) = (2,)
function Base.getindex(E::Edge, i::Int)
  if i == 1
    return E.p1
  elseif i == 2
    return E.p2
  else
    BoundsError(E, i)
  end
end)
