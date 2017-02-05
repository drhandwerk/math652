__precompile__()
module UniformMesh

using PyPlot

export UniformTriangleMesh, drawmesh
"""
  UniformTriangleMesh(m,n)

Creates a uniform triangle mesh of the unit square with 'm' partitions in the
x direction and 'n' partitions in the y direction.
"""
type UniformTriangleMesh
  m::Int64 # Number of partitions in x-dir
  n::Int64 # Number of partitions in y-dir
  num_vertices::Int64
  vertices::Array{Float64, 2}
  triangles::Array{Float64, 2}
  edges::Array{Float64, 2}

  function UniformTriangleMesh(m::Int64,n::Int64)
    mesh = new(m,
               n,
               (m+1)*(n+1),
               Array{Float64, 2}((m+1)*(n+1), 2),
               Array{Float64, 2}(m*n*2, 6),
               Array{Float64, 2}(3*m*n + n + m, 4))
    generateVertices!(mesh)
    generateTriangles!(mesh)
    generateEdges!(mesh)
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
      mesh.vertices[j + (i-1)*(mesh.n+1), 1] = width*(i-1)  # x component
      mesh.vertices[j + (i-1)*(mesh.n+1), 2]=  height*(j-1) # y component
    end
  end
end

"""
  generateTriangles!(mesh)

Creates all of the triangles for the mesh overwriting the Constructor created
empty array in place. Numbering starts at bottom left, proceeds up, then right.
"""
function generateTriangles!(mesh::UniformTriangleMesh)
  count::Int32 = 1
  for i = 1:mesh.m
    for j = 1:mesh.n
      base::Int64 = j + (mesh.n + 1)*(i-1)
      mesh.triangles[count, 1:2] = mesh.vertices[base, :]
      mesh.triangles[count, 3:4] = mesh.vertices[base + (mesh.n + 1), :]
      mesh.triangles[count, 5:6] = mesh.vertices[base + (mesh.n + 2), :]
      count += 1
      mesh.triangles[count, 1:2] = mesh.vertices[base, :]
      mesh.triangles[count, 3:4] = mesh.vertices[base + (mesh.n + 2), :]
      mesh.triangles[count, 5:6] = mesh.vertices[base + 1, :]
      count += 1
    end
  end
end

"""
  generateEdges!(mesh)

Creates all of the edges for the mesh, overwriting the default Constructor
generated Array in the process. First count vertical edges, then horizontal edges,
then diagonal edges.
"""
function generateEdges!(mesh::UniformTriangleMesh)
  count::Int32 = 1
  # Vertical Edges
  for j = 1:mesh.m+1
    for i = 1:mesh.n
      mesh.edges[count, 1:2] = mesh.vertices[i + (mesh.n + 1)*(j - 1), :]
      mesh.edges[count, 3:4] = mesh.vertices[(i + 1) + (mesh.n + 1)*(j - 1), :]
      count += 1
    end
  end
  # Horizontal Edges
  for j = 1:mesh.m
    for i = 1:mesh.n+1
      mesh.edges[count, 1:2] = mesh.vertices[i + (mesh.n + 1)*(j - 1), :]
      mesh.edges[count, 3:4] = mesh.vertices[(i + mesh.n + 1) + (mesh.n + 1)*(j - 1), :]
      count += 1
    end
  end
  # Diagonal Edges
  for j = 1:mesh.m
    for i = 1:mesh.n
      mesh.edges[count, 1:2] = mesh.vertices[i + (mesh.n + 1)*(j - 1), :]
      mesh.edges[count, 3:4] = mesh.vertices[(i + mesh.n + 2) + (mesh.n + 1)*(j - 1), :]
      count += 1
    end
  end
end

"""
  drawmesh(mesh)

Plot the vertices and edges for the mesh.
"""
function drawmesh(mesh::UniformTriangleMesh)
  # Vertices
  scatter(mesh.vertices[:,1], mesh.vertices[:,2], marker="o", s = 100, color="blue")
  # Edges
  plot([mesh.edges[:,1],mesh.edges[:,3]],[mesh.edges[:,2],mesh.edges[:,4]],linestyle="-",color="green")
end

end # End module
