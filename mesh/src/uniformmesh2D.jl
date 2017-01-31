# Build uniform triangluar mesh over unit square
include("vertex.jl")
include("triangle.jl")

type Uniformmesh2D
  m::Int64 # Number of partitions in x-dir
  n::Int64 # Number of partitions in y-dir
  num_vertices::Int64
  verticies::Array{Vertex}
  triangles::Array{Triangle}

  function Uniformmesh2D(m::Int64,n::Int64)
    this = new(m,n,(m+1)*(n+1),Array{Vertex}((m+1)*(n+1)),Array{Triangle}(m*n*2))
    generateVertices(this)
    generateTriangles(this)
    return this
  end
end

# Generate array of verticies numbered from bottom left. Increaes up then right.
function generateVertices(mesh::Uniformmesh2D)
  width::Float64 = 1.0/mesh.m
  height::Float64 = 1.0/mesh.n
  for i = 1:mesh.m+1
    for j = 1:mesh.n+1
      mesh.verticies[j + (i-1)*(mesh.n+1)] = Vertex(width*(i-1),height*(j-1))
      println(i,j)
    end
  end
end
# Generate array of triangles numbered from bottom left. Increases up then right.
function generateTriangles(mesh::Uniformmesh2D)

end
