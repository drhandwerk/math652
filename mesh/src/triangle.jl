include("Vertex.jl")

"""
  Triangle(p1, p2, p3)

Creates a Triangle with vertices 'p1', 'p2', 'p3'.
"""
type Triangle <: AbstractArray{Float64, 1}
  p1::Vertex
  p2::Vertex
  p3::Vertex
  function Triangle(p1::Vertex,p2::Vertex,p3::Vertex)
    triangle = new(p1,p2,p3)
    #sort!(triangle) #TODO if not using sort! can put on one line
    return triangle
  end
end
Base.size(T::Triangle) = (3,)
function Base.getindex(T::Triangle, i::Int)
  if i == 1
    return T.p1
  elseif i == 2
    return T.p2
  elseif i == 3
    return T.p3
  else
    BoundsError(T, i)
  end
end

"""
  sort!(T)

Arranges the verticies of the triangle in a counter-clockwise order. First find
the center of the triangle and then compute the angle between the center and all
the other verticies. Arrange in increasing order.

Numbering starts with the bottom left vertex.
"""
function Base.sort!(T::Triangle)
  temp_triangle::Triangle = deepcopy(T)
  center::Vector = (T.p1 + T.p2 + T.p3)./3
  angles::Array{Float64} = [atan(dot(center, T.p1)),
                            atan(dot(center, T.p2)),
                            atan(dot(center, T.p3))]
  p1_index::Int = indmin(angles)
  p3_index::Int = indmax(angles)
  p2_index::Array{Int} = setdiff([1,2,3],[p1_index, p3_index])
  T.p1 = temp_triangle[p1_index]
  T.p2 = temp_triangle[p2_index[1]]
  T.p3 = temp_triangle[p3_index]
end
