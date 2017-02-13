__precompile__()
module UniformMesh

using PyPlot

export UniformTriangleMesh
export UniformRectMesh
export drawmesh

abstract UniformPolyMesh

type UniformRectMesh <: UniformPolyMesh
  m::Int64
  n::Int64
  vertices::Array{Float64,2}
  rectangles::Array{Int64,2}
  edges::Array{Int64,2}

  function UniformRectMesh(m::Int64,n::Int64)
    mesh = new(m,
               n,
               Array{Float64, 2}((m+1)*(n+1), 2),
               Array{Int64, 2}(m*n, 4),
               Array{Int64, 2}(2*m*n + n + m, 2))
    generateVertices!(mesh)
    generateRectangles!(mesh)
    generateEdges!(mesh)
    return mesh
  end
end

"""
  UniformTriangleMesh(m,n)

Creates a uniform triangle mesh of the unit square with 'm' partitions in the
x direction and 'n' partitions in the y direction.

TODO: should triangles and edges be points (currently) or the indices of the
vertices? e.g.
mesh.triangles[count, 1:2] = mesh.vertices[base, :]
mesh.triangles[count, 3:4] = mesh.vertices[base + (mesh.n + 1), :]
mesh.triangles[count, 5:6] = mesh.vertices[base + (mesh.n + 2), :]

or

mesh.triangles[count, 1] = base
mesh.triangles[count, 2] = base + (mesh.n + 1)
mesh.triangles[count, 3] = base + (mesh.n + 2)
Don't forget to change constructor if this is implemented
"""
type UniformTriangleMesh <: UniformPolyMesh
  m::Int64 # Number of partitions in x-dir
  n::Int64 # Number of partitions in y-dir
  vertices::Array{Float64, 2}
  triangles::Array{Int64, 2}
  edges::Array{Int64, 2}

  function UniformTriangleMesh(m::Int64,n::Int64)
    mesh = new(m,
               n,
               Array{Float64, 2}((m+1)*(n+1), 2),
               Array{Int64, 2}(m*n*2, 3),
               Array{Int64, 2}(3*m*n + n + m, 2))
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
function generateVertices!(mesh::UniformPolyMesh)
  width::Float64 = 1.0/mesh.m
  height::Float64 = 1.0/mesh.n
  for i = 1:mesh.m+1
    for j = 1:mesh.n+1
      mesh.vertices[j + (i-1)*(mesh.n+1), 1] = width*(i-1)  # x component
      mesh.vertices[j + (i-1)*(mesh.n+1), 2] = height*(j-1) # y component
    end
  end
end

"""
  generateRectangles!(UniformRectMesh)

Creates all of the rectangles(indices of vertices) for the mesh overwriting the Constructor created
empty array in place. Numbering starts at bottom left, proceeds up, then right.
"""
function generateRectangles!(mesh::UniformRectMesh)
  count::Int32 = 1
  for i = 1:mesh.m
    for j = 1:mesh.n
      base::Int64 = j + (mesh.n + 1)*(i-1)
      mesh.rectangles[count, 1] = base
      mesh.rectangles[count, 2] = base + (mesh.n + 1)
      mesh.rectangles[count, 3] = base + (mesh.n + 2)
      mesh.rectangles[count, 4] = base + 1
      count += 1
    end
  end
end

"""
  generateTriangles!(mesh)

Creates all of the triangles(indices of vertices) for the mesh overwriting the Constructor created
empty array in place. Numbering starts at bottom left, proceeds up, then right.
"""
function generateTriangles!(mesh::UniformTriangleMesh)
  count::Int32 = 1
  for i = 1:mesh.m
    for j = 1:mesh.n
      base::Int64 = j + (mesh.n + 1)*(i-1)
      # two triangles to maintain vertex ordering
      mesh.triangles[count, 1] = base
      mesh.triangles[count, 2] = base + (mesh.n + 1)
      mesh.triangles[count, 3] = base + (mesh.n + 2)
      count += 1
      mesh.triangles[count, 1] = base
      mesh.triangles[count, 2] = base + (mesh.n + 2)
      mesh.triangles[count, 3] = base + 1
      count += 1
    end
  end
end

"""
  generateEdges!(UniformRectMesh)

Creates all of the edges for the mesh, overwriting the default Constructor
generated Array in the process. First count vertical edges, then horizontal edges.
"""
function generateEdges!(mesh::UniformRectMesh)
  count::Int32 = 1
  # Vertical Edges
  for j = 1:mesh.m+1
    for i = 1:mesh.n
      mesh.edges[count, 1] = i + (mesh.n + 1)*(j - 1)
      mesh.edges[count, 2] = (i + 1) + (mesh.n + 1)*(j - 1)
      count += 1
    end
  end
  # Horizontal Edges
  for j = 1:mesh.m
    for i = 1:mesh.n+1
      mesh.edges[count, 1] = i + (mesh.n + 1)*(j - 1)
      mesh.edges[count, 2] = (i + mesh.n + 1) + (mesh.n + 1)*(j - 1)
      count += 1
    end
  end
end

"""
  generateEdges!(UniformTriangleMesh)

Creates all of the edges for the mesh, overwriting the default Constructor
generated Array in the process. First count vertical edges, then horizontal edges,
then diagonal edges.
"""
function generateEdges!(mesh::UniformTriangleMesh)
  count::Int32 = 1
  # Vertical Edges
  for j = 1:mesh.m+1
    for i = 1:mesh.n
      mesh.edges[count, 1] = i + (mesh.n + 1)*(j - 1)
      mesh.edges[count, 2] = (i + 1) + (mesh.n + 1)*(j - 1)
      count += 1
    end
  end
  # Horizontal Edges
  for j = 1:mesh.m
    for i = 1:mesh.n+1
      mesh.edges[count, 1] = i + (mesh.n + 1)*(j - 1)
      mesh.edges[count, 2] = (i + mesh.n + 1) + (mesh.n + 1)*(j - 1)
      count += 1
    end
  end
  # Diagonal Edges
  for j = 1:mesh.m
    for i = 1:mesh.n
      mesh.edges[count, 1] = i + (mesh.n + 1)*(j - 1)
      mesh.edges[count, 2] = (i + mesh.n + 2) + (mesh.n + 1)*(j - 1)
      count += 1
    end
  end
end

"""
  drawmesh(mesh)

Plot the vertices and edges for the mesh.
"""
function drawmesh(mesh::UniformPolyMesh)
  if mesh.n > 32 || mesh.m > 32
    throw(ArgumentError("Drawing a fine mesh(>32,>32) is a bad idea."))
  end
  # Edges
  for i = 1:size(mesh.edges,1)
  plot(mesh.vertices[mesh.edges[i,:],1], mesh.vertices[mesh.edges[i,:],2],linestyle="-",color="green")
  end
  # Vertices
  scatter(mesh.vertices[:,1], mesh.vertices[:,2], marker="o", s = 30, color="blue")
  axis("square")
end

end # End module
