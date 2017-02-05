module UniformMesh

"""
  UniformMesh(n,m)

Abstract type for a uniform mesh. Subtypes could include triangle or rectangle.
"""
abstract UniformMesh
  m::Int64 # Number of partitions in x-dir
  n::Int64 # Number of partitions in y-dir
  num_vertices::Int64
  vertices::Array{Float64, 2}
  triangles::Array{Float64, 2}
  edges::Array{Float64, 2}
end

end # End module
