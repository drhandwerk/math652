# Build uniform triangluar mesh over unit square
include("vertex.jl")
include("triangle.jl")

type Uniformmesh2D
  m::Int64 # Number of partitions in x-dir
  n::Int64 # Number of partitions in y-dir
  num_vertices::Int64
  verticies::Array{Vertex}
  triangles::Array{Triangle}

  function Uniformmesh2D(m,n)
    this = new(m,n,(m+1)*(n+1))
    generateVertices(this.mesh)
    return this
  end
end

# Generate array of verticies numbered from bottom left. Increaes up then right.
function generateVertices(mesh::Uniformmesh2D)
  println(mesh.n)
  mesh.verticies = Array{Triangle}(mesh.m*mesh.n)

  for i = 1:mesh.n+1
    for j = 0:mesh.m
      mesh.verticies[j + (i-1)*(mesh.m+1)] = Vertex(0.,0.)
    end
  end
end
# Generate array of triangles numbered from bottom left. Increases up then right.
function generateTriangles()

end
