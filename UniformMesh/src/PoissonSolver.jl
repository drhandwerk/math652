
include("UniformMesh.jl")
using UniformMesh
using Cubature # for quadrature

export poissonsolve

function poissonsolve(n::Int64)
  # get mesh
  mesh = UniformRectMesh(n,n)
  # get global stiffness matrix
  G = gsm(mesh)
  # set RHS from quadrature over elements
  b = rhs(mesh)
  # set BC
  setdirichlet!(mesh,G,b)
  # compute coeffs
  c = G\b
  # return as matrix
  return reshape(c,n+1,n+1)
end

"""
  elementrhs()
Do quadrature over an element for RHS.
"""
function elementrhs(length::Int64, p1::Array{Float64,1}, p2::Array{Float64,1})
  area = (p2[1] - p1[1]) * (p2[2] - p1[2])
  (b1, err1) = hcubature(x -> 2*pi^2*sin(pi*x[1])*sin(pi*x[2]) * ((p2[1] - x[1])*(p2[2] - x[2]))/area, p1, p2)
  (b2, err2) = hcubature(x -> 2*pi^2*sin(pi*x[1])*sin(pi*x[2]) * ((x[1] - p1[1])*(p2[2] - x[2]))/area, p1, p2)
  (b3, err3) = hcubature(x -> 2*pi^2*sin(pi*x[1])*sin(pi*x[2]) * ((x[1] - p1[1])*(x[2] - p1[2]))/area, p1, p2)
  (b4, err4) = hcubature(x -> 2*pi^2*sin(pi*x[1])*sin(pi*x[2]) * ((p2[1] - x[1])*(x[2] - p1[2]))/area, p1, p2)
  return [b1;b2;b3;b4]
end

"""
  rhs(mesh)
Assemble the vector b (RHS) by summing over elements.
"""
function rhs(mesh::UniformRectMesh)
  b = zeros(Float64, size(mesh.vertices,1), 1) # preallocate
  for i = 1:size(mesh.rectangles,1)
    @inbounds elemb = elementrhs(size(mesh.vertices,1),
                       mesh.vertices[mesh.rectangles[i,1],:],
                       mesh.vertices[mesh.rectangles[i,3],:])
    b[mesh.rectangles[i,:]] += elemb
  end
  return b
end
"""
  esmrect(p1, p2)
Send in lower left corner 'p1' and upper right corner 'p2' of a rectangle. Return the
element stiffness matrix.
"""
function esmrect(p1::Array{Float64,1}, p2::Array{Float64,1})
  dx = p2[1] - p1[1]
  dy = p2[2] - p1[2]
  E = zeros(Float64, 4, 4)
  A11 = 1./3. * (dx/dy + dy/dx)
  A12 = 1./6. * dx/dy - 1./3. * dy/dx
  A13 = -1./6. * (dx/dy + dy/dx)
  A14 = -1./3. * dx/dy + 1./6. * dy/dx

  E = [A11 A12 A13 A14;
       A12 A11 A14 A13;
       A13 A14 A11 A12;
       A14 A13 A12 A11]

end

"""
  gsm(mesh)
Assemble the gsm from all of the element stiffness matrices.
"""
function gsm(mesh::UniformRectMesh)
  G = zeros(Float64, size(mesh.vertices,1), size(mesh.vertices,1)) # preallocate
  #G = SharedArray(Float64, (size(mesh.vertices,1), size(mesh.vertices,1)))
  for i = 1:size(mesh.rectangles,1)
    E = esmrect(mesh.vertices[mesh.rectangles[i,1],:],
                mesh.vertices[mesh.rectangles[i,3],:])
    G[mesh.rectangles[i,:], mesh.rectangles[i,:]] += E
  end
  return G
end


"""
  setdirichlet!(mesh, G, b)
Naively modifies in place the esm, G and the RHS, b for homogeneous Dirichlet BC.
"""
function setdirichlet!(mesh::UniformRectMesh, G::Array{Float64, 2}, b::Array{Float64,2})
  #b = zeros(size(mesh.vertices,1), 1) # preallocate
  # get boundary vertices
  V = reshape(1:size(mesh.vertices,1), mesh.m+1, mesh.n+1)
  interiorvertices = V[2:end-1, 2:end-1][:] # vector
  exteriorvertices = setdiff(1:size(mesh.vertices,1), interiorvertices)

  ind = trues(size(G)) # indices to zero out
  ind[interiorvertices,:] = false # keep interiorvertices rows
  for i in exteriorvertices # keep what's on the diagonal
    @inbounds ind[i,i] = false
  end
  @inbounds G[ind] = 0.0
  #= this is to set the Dirichlet boundary to some function.
  X,Y = ndgrid(linspace(0,1,mesh.n+1),linspace(0,1,mesh.m+1))
  g = X + Y
  bvals = g[exteriorvertices];
  =#
  @inbounds b[exteriorvertices] = 0.0
end

function ndgrid{T}(v1::AbstractArray{T}, v2::AbstractArray{T})
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repmat(v1, 1, n), repmat(v2, m, 1))
end
