# Build uniform triangluar mesh over unit square

#TODO use modules
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
  triangles::Array{Array{Float64, 2}}
  edges::Array{Array{Float64, 2}}

  function UniformTriangleMesh(m::Int64,n::Int64)
    mesh = new(m,
               n,
               (m+1)*(n+1),
               Array{Float64, 2}((m+1)*(n+1), 2),
               Array{Array{Float64, 2}}(m*n*2),
               Array{Array{Float64, 2}}(3*m*n + n + m))
    generateVertices!(mesh)
    #generateTriangles!(mesh)
    #generateEdges!(mesh)
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
      mesh.vertices[j + (i-1)*(mesh.n+1), 1] = width*(i-1)
      mesh.vertices[j + (i-1)*(mesh.n+1), 2]=  height*(j-1)
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
      mesh.triangles[count] = [mesh.vertices[base]
                               mesh.vertices[base + (mesh.n + 1)]
                               mesh.vertices[base + (mesh.n + 2)]]
      count += 1
      mesh.triangles[count] = [mesh.vertices[base]
                               mesh.vertices[base + (mesh.n + 2)]
                               mesh.vertices[base + 1]]
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
      mesh.edges[count] = [mesh.vertices[i + (mesh.n + 1)*(j - 1)]
                           mesh.vertices[(i + 1) + (mesh.n + 1)*(j - 1)]]
      count += 1
    end
  end
  # Horizontal Edges
  for j = 1:mesh.m
    for i = 1:mesh.n+1
      mesh.edges[count] = [mesh.vertices[i + (mesh.n + 1)*(j - 1)]
                           mesh.vertices[(i + mesh.n + 1) + (mesh.n + 1)*(j - 1)]]
      count += 1
    end
  end
  # Diagonal Edges
  for j = 1:mesh.m
    for i = 1:mesh.n
      mesh.edges[count] = [mesh.vertices[i + (mesh.n + 1)*(j - 1)]
                           mesh.vertices[(i + mesh.n + 2) + (mesh.n + 1)*(j - 1)]]
      count += 1
    end
  end
end
