__precompile__()
module PoissonSolverTri
@everywhere push!(LOAD_PATH,"./")
include("UniformMesh.jl")

using UniformMesh
using Cubature # for quadrature
using PyPlot

export poissonsolve
export solveanddraw

function poissonsolve(n::Int64)
  # get mesh
  mesh = UniformTriangleMesh(n,n)
  # get global stiffness matrix
  G = gsm(mesh)
  # set RHS from quadrature over elements
  b = rhs(mesh)
  # set BC
  setneumann!(mesh,G,b)
  setdirichlet!(mesh,G,b)
  # compute coeffs
  c = G\b
  # return as matrix
  return reshape(c,n+1,n+1)
end


"""
  solveanddraw(n)
Solves the problem and plots the solution.
"""
function solveanddraw(n::Int64)
  c = poissonsolve(n)
  trimesh = UniformTriangleMesh(n,n)
  fig = figure()
  ax = fig[:add_subplot](111, projection="3d")
  ax[:plot_trisurf](trimesh.vertices[:,1],trimesh.vertices[:,2], triangles=trimesh.triangles-1, c[:], alpha=1, cmap="viridis")
  ax[:set_title](string(n, "x", n), fontsize=22)
  ax[:set_xlabel]("X", fontsize=20)
  ax[:set_ylabel]("Y", fontsize=20)
  ax[:xaxis][:set_tick_params](labelsize=18)
  ax[:yaxis][:set_tick_params](labelsize=18)
  ax[:zaxis][:set_tick_params](labelsize=18)
end

"""
  elementrhs()
Do quadrature over a triangular element for RHS.
"""
function elementrhs(length::Int64, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  b1 = trigaussquad(x -> 2*pi^2*sin(pi*x[1])*sin(pi*x[2]) * 0.5 * abs(det([x[1] x[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))/area, p1, p2, p3)
  b2 = trigaussquad(x -> 2*pi^2*sin(pi*x[1])*sin(pi*x[2]) * 0.5 * abs(det([p1[1] p1[2] 1; x[1] x[2] 1; p3[1] p3[2] 1]))/area, p1, p2, p3)
  b3 = trigaussquad(x -> 2*pi^2*sin(pi*x[1])*sin(pi*x[2]) * 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; x[1] x[2] 1]))/area, p1, p2, p3)
  return [b1; b2; b3]
end

"""
  trigaussquad()
Do Gaussian quadrature over a triangle using three interpolation points.
"""
function trigaussquad(f::Function, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  barycoords = [2.0/3.0 1.0/6.0 1.0/6.0;
                1.0/6.0 2.0/3.0 1.0/6.0;
                1.0/6.0 1.0/6.0 2.0/3.0]
  sum = 0.0
  for i = 1:3
    sum += 1.0/3.0 * f(barycoords[i,1]*p1 + barycoords[i,2]*p2 + barycoords[i,3]*p3)
  end
  sum *= area
  return sum
end
"""
  rhs(mesh)
Assemble the vector b (RHS) by summing over elements.
"""
function rhs(mesh::UniformTriangleMesh)
  b = zeros(Float64, size(mesh.vertices,1), 1) # preallocate
  for i = 1:size(mesh.triangles,1)
    @inbounds elemb = elementrhs(size(mesh.vertices,1),
                       mesh.vertices[mesh.triangles[i,1],:],
                       mesh.vertices[mesh.triangles[i,2],:],
                       mesh.vertices[mesh.triangles[i,3],:])
    b[mesh.triangles[i,:]] += elemb
  end
  return b
end
"""
  esmtri(p1, p2)
Send in three vertices 'p1', 'p2', 'p3' of a triangle. Return the
element stiffness matrix.
"""
function esmtri(p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  E = zeros(Float64, 3, 3)
  A11 = 1.0/(2.0*area) * (p2[2] - p3[2])^2 + (p1[1] - p3[1])^2
  A22 = 1.0/(2.0*area) * (p3[2] - p1[2])^2 + (p1[1] - p3[1])^2
  A33 = 1.0/(2.0*area) * (p1[2] - p2[2])^2 + (p2[1] - p1[1])^2
  A12 = 1.0/(2.0*area) * (p2[2] - p3[2])*(p3[2] - p1[2]) + (p1[1] - p3[1])*(p1[1] - p3[1])
  A13 = 1.0/(2.0*area) * (p2[2] - p3[2])*(p1[2] - p2[2]) + (p1[1] - p3[1])*(p2[1] - p1[1])
  A23 = 1.0/(2.0*area) * (p3[2] - p1[2])*(p1[2] - p2[2]) + (p1[1] - p3[1])*(p2[1] - p1[1])

  E = [A11 A12 A13;
       A12 A22 A23;
       A13 A23 A33]

end

"""
  gsm(mesh)
Assemble the gsm from all of the element stiffness matrices.
"""
function gsm(mesh::UniformTriangleMesh)
  G = zeros(Float64, size(mesh.vertices,1), size(mesh.vertices,1)) # preallocate
  #G = SharedArray(Float64, (size(mesh.vertices,1), size(mesh.vertices,1)))
  for i = 1:size(mesh.triangles,1)
    E = esmtri(mesh.vertices[mesh.triangles[i,1],:],
                mesh.vertices[mesh.triangles[i,2],:],
                mesh.vertices[mesh.triangles[i,3],:])
    G[mesh.triangles[i,:], mesh.triangles[i,:]] += E
  end
  return G
end


"""
  setneumann!(mesh, G, b)
Modify b using quadrature on the boundary.
"""
function setneumann!(mesh::UniformTriangleMesh, G::Array{Float64, 2}, b::Array{Float64, 2})
  # get top and bottom elements
  topelements = zeros(Int64, mesh.m)
  bottomelements = zeros(Int64, mesh.m)
  count = 1
  for i = 1:mesh.m
    topelements[count] = i*2*mesh.n
    count += 1
  end
  count = 1
  for i = 0:mesh.m-1
    bottomelements[count] = i*2*mesh.n + 1
    count += 1
  end
  # quadrature over each elements
  for e in topelements
    p1 = mesh.vertices[mesh.triangles[e,1],:]
    p2 = mesh.vertices[mesh.triangles[e,2],:]
    p3 = mesh.vertices[mesh.triangles[e,3],:]
    area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
    (b1, err1) = hcubature(x -> pi*sin(pi*x[1])*cos(pi*x[2]) * 0.5 * abs(det([p1[1] p1[2] 1; x[1] x[2] 1; p3[1] p3[2] 1]))/area, p3, p2)
    (b2, err1) = hcubature(x -> pi*sin(pi*x[1])*cos(pi*x[2]) * 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; x[1] x[2] 1]))/area, p3, p2)
    b[mesh.triangles[e,:]] += [0.0; b1; b2]
  end
  for e in bottomelements
    p1 = mesh.vertices[mesh.triangles[e,1],:]
    p2 = mesh.vertices[mesh.triangles[e,2],:]
    p3 = mesh.vertices[mesh.triangles[e,3],:]
    area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
    (b1, err1) = hcubature(x -> -pi*sin(pi*x[1])*cos(pi*x[2]) * 0.5 * abs(det([x[1] x[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))/area, p1, p2)
    (b2, err1) = hcubature(x -> -pi*sin(pi*x[1])*cos(pi*x[2]) * 0.5 * abs(det([p1[1] p1[2] 1; x[1] x[2] 1; p3[1] p3[2] 1]))/area, p1, p2)
    b[mesh.triangles[e,:]] += [b1; b2; 0.0]
  end
  # update global b with the new parts from quadrature
end

"""
  setdirichlet!(mesh, G, b)
Naively modifies in place the esm, G and the RHS, b for homogeneous Dirichlet BC on left and right.
"""
function setdirichlet!(mesh::UniformTriangleMesh, G::Array{Float64, 2}, b::Array{Float64,2})
  # get all boundary vertices
  V = reshape(1:size(mesh.vertices,1), mesh.m+1, mesh.n+1)
  interiorvertices = V[2:end-1, 2:end-1][:] # vector
  exteriorvertices = setdiff(1:size(mesh.vertices,1), interiorvertices)
  topvertices = V[1,2:end-1][:] # only need middle of top
  bottomvertices = V[end,2:end-1][:] # only need middle of bottom
  leftvertices = V[:,1][:]
  rightvertices = V[:,end][:]

  ind = trues(size(G)) # indices to zero out
  ind[interiorvertices,:] = false # keep interiorvertices rows
  ind[topvertices,:] = false # keep topvertices
  ind[bottomvertices,:] = false # keep bottomvertices
  for i in leftvertices # keep what's on the diagonal
    @inbounds ind[i,i] = false
  end
  for i in rightvertices
    @inbounds ind[i,i] = false
  end
  @inbounds G[ind] = 0.0
  #= this is to set the Dirichlet boundary to some function.
  X,Y = ndgrid(linspace(0,1,mesh.n+1),linspace(0,1,mesh.m+1))
  g = X + Y
  bvals = g[exteriorvertices];
  =#
  @inbounds b[leftvertices] = 0.0
  @inbounds b[rightvertices] = 0.0
end

function ndgrid{T}(v1::AbstractArray{T}, v2::AbstractArray{T})
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repmat(v1, 1, n), repmat(v2, m, 1))
end

end # module
