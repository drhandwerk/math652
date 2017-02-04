"""
  Vertex(x, y)

Creates a two-dimensional vertex point ('x','y'). Inherits AbstractArray for all the built in Vector
functions. Entries are floats.

"""
type Vertex <: AbstractArray{Float64, 1}
  x::Float64
  y::Float64
end
# These functions are needed for the inheritence.
Base.size(V::Vertex) = (2,)
function Base.getindex(V::Vertex, i::Int)
  if i == 1
    return V.x
  elseif i == 2
    return V.y
  else
    BoundsError(V, i)
  end
end
