# Build uniform triangluar mesh over unit square

#TODO use modules

include("Vertex.jl")
include("Triangle.jl")
"""
  UniformTriangleMesh(m,n)

Creates a uniform triangle mesh of the unit square with 'm' partitions in the
x direction and 'n' partitions in the y direction.
"""
type UniformTriangleMesh
  m::Int64 # Number of partitions in x-dir
  n::Int64 # Number of partitions in y-dir
  num_vertices::Int64
  vertices::Array{Vertex}
  triangles::Array{Triangle}

  function UniformTriangleMesh(m::Int64,n::Int64)
    mesh = new(m,n,(m+1)*(n+1),Array{Vertex}((m+1)*(n+1)),Array{Triangle}(m*n*2))
    generateVertices!(mesh)
    generateTriangles!(mesh)
    return mesh
  end
end

"""
  generateVertices!(mesh)

Creates all of the vertices for the mesh and changes in place. Numbering starts
at the bottom left, and goes up each column before going to the bottom of the next
column.
"""
function generateVertices!(mesh::UniformTriangleMesh)
  width::Float64 = 1.0/mesh.m
  height::Float64 = 1.0/mesh.n
  for i = 1:mesh.m+1
    for j = 1:mesh.n+1
      mesh.vertices[j + (i-1)*(mesh.n+1)] = Vertex(width*(i-1),height*(j-1))
    end
  end
end

"""
  generateTriangles!(mesh)

Creates all of the triangles for the mesh overwriting the Constructor created
empty array in place. Numbering starts at bottom left, proceeds up, then right.
"""
function generateTriangles!(mesh::UniformTriangleMesh)
  for i = 1:mesh.m
    for j = 1:2*mesh.n
      #TODO fix which vertices get used for triangles
      if j % 2 == 1
        mesh.triangles[j + (i-1)*mesh.n*2] = Triangle(mesh.vertices[i],
                                                      mesh.vertices[i+(mesh.n+1)],
                                                      mesh.vertices[i+(mesh.n+2)])
      else
        mesh.triangles[j + (i-1)*mesh.n*2] = Triangle(mesh.vertices[i],
                                                      mesh.vertices[i+(mesh.n+2)],
                                                      mesh.vertices[i+1])
      end
    end
  end
end
