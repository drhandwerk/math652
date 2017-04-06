__precompile__()
module RitzProjection
@everywhere push!(LOAD_PATH,"./")
include("UniformMesh.jl")

using UniformMesh
using PyPlot

"""
  function ritzerrors()
Calculate the ritz errors for finer and finer meshes.
"""
function ritzerrors()
  u(x) = sin(pi*x[1]).*sin(pi*x[2])
  Du(x) = [pi*cos(pi*x[1])*sin(pi*x[2]), pi*sin(pi*x[1])*cos(pi*x[2])]
  println(" ", "="^43)
  println("|", " "^5, "n", " "^3, "|", " "^4, "L2error", " "^5, "|",
          " "^5, "L2Conv", " "^5, "|")
  println(" ", "="^43)
  prevL2error = 1.0
  for n = [8,16,32,64]
    mesh = UniformTriangleMesh(n,n)
    errL2 = ritzL2error(mesh, u, Du)
    convergenceL2 = log(2, prevL2error/errL2)
    prevL2error = errL2
    print("|", " "^3)
    @printf("%3d", n)
    print(" "^3, "|", " "^3)
    @printf("%6.8f", errL2)
    print(" "^3, "|", " "^3)
    @printf("%6.8f", convergenceL2)
    println(" "^3, "|", " "^3)
  end
  println(" ", "="^43)
end

"""
  ritzL2error(mesh, f)
Calculate the L2 norm error of the ritz projection.
"""
function ritzL2error(mesh::UniformTriangleMesh, f::Function, Df::Function)
  B = gsm(mesh) # global stiffness matrix
  b = ritzrhs(mesh, Df) # load vector
  setalldirichlet!(mesh, B, b)
  c = B\b # ritz projection solution
  #fig = figure()
  #ax = fig[:add_subplot](111, projection="3d")
  #ax[:plot_trisurf](mesh.vertices[:,1],mesh.vertices[:,2], triangles=mesh.triangles-1, c[:], alpha=.8, cmap="viridis", edgecolors=:black)
  errL2 = errorL2(c, mesh.n, f)
end

"""
  elementrhs()
Do quadrature over a triangular element for RHS.
"""
function elementritzrhs(p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1}, Df::Function)
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  b1 = trigaussquad(x -> vecdot(Df(x), [(p2[2] - p3[2]), (p3[1] - p2[1])]./(2*area)), p1, p2, p3)
  b2 = trigaussquad(x -> vecdot(Df(x), [(p3[2] - p1[2]), (p1[1] - p3[1])]./(2*area)), p1, p2, p3)
  b3 = trigaussquad(x -> vecdot(Df(x), [(p1[2] - p2[2]), (p2[1] - p1[1])]./(2*area)), p1, p2, p3)
  return [b1; b2; b3]
end

"""
  rhs(mesh)
Assemble the vector b (RHS) by summing over elements.
"""
function ritzrhs(mesh::UniformTriangleMesh, Df::Function)
  b = zeros(Float64, size(mesh.vertices,1), 1) # preallocate
  for i = 1:size(mesh.triangles,1)
    elemb = elementritzrhs(mesh.vertices[mesh.triangles[i,1],:],
                       mesh.vertices[mesh.triangles[i,2],:],
                       mesh.vertices[mesh.triangles[i,3],:],
                       Df)
    b[mesh.triangles[i,:]] += elemb
  end
  return b
end


"""
  errorL2(c, n, u)
Compute the L2 error between the exact solution 'u' and the approximate solution 'c'.
"""
function errorL2(c::Array{Float64, 2}, n::Int64, u::Function)
  mesh = UniformTriangleMesh(n,n)
  error = 0.0
  for t = 1:size(mesh.triangles,1)
    trierror = 0.0
    p1 = mesh.vertices[mesh.triangles[t,1],:]
    p2 = mesh.vertices[mesh.triangles[t,2],:]
    p3 = mesh.vertices[mesh.triangles[t,3],:]
    area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
    c1,c2,c3 = c[:][mesh.triangles[t,:]]
    uh(x) = x == p1 ? c1 : x == p2 ? c2 : c3
    # integrate the difference of the exact and approximate solutions
    trierror = trigaussquad(x -> norm(u(x) - uh(p1) * 0.5* abs(det([x[1] x[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))/area
                                           - uh(p2) * 0.5* abs(det([p1[1] p1[2] 1; x[1] x[2] 1; p3[1] p3[2] 1]))/area
                                           - uh(p3) * 0.5* abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; x[1] x[2] 1]))/area)^2, p1, p2, p3)

    error += trierror
  end
  return sqrt(error)
end

"""
  esmtri(p1, p2)
Send in three vertices 'p1', 'p2', 'p3' of a triangle. Return the
element stiffness matrix.
"""
function esmtri(p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  E = zeros(Float64, 3, 3)
  A11 = 1.0/(4.0*area) * ((p2[2] - p3[2])^2 + (p3[1] - p2[1])^2)
  A22 = 1.0/(4.0*area) * ((p3[2] - p1[2])^2 + (p1[1] - p3[1])^2)
  A33 = 1.0/(4.0*area) * ((p1[2] - p2[2])^2 + (p2[1] - p1[1])^2)
  A12 = 1.0/(4.0*area) * ((p2[2] - p3[2])*(p3[2] - p1[2]) + (p3[1] - p2[1])*(p1[1] - p3[1]))
  A13 = 1.0/(4.0*area) * ((p2[2] - p3[2])*(p1[2] - p2[2]) + (p3[1] - p2[1])*(p2[1] - p1[1]))
  A23 = 1.0/(4.0*area) * ((p3[2] - p1[2])*(p1[2] - p2[2]) + (p1[1] - p3[1])*(p2[1] - p1[1]))

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
  trigaussquad(f,p1,p2,p3)
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
  return sum * area
end

"""
  lingaussquad(f,p1,p2)
Perform Gaussian quadrature over a line using three interpolation points.
"""
function lingaussquad(f::Function, p1::Array{Float64,1}, p2::Array{Float64,1})
  length = sqrt((p2[1]-p1[1])^2 + (p2[2]-p1[2])^2)
  weights = [5./18. 8./18. 5./18.]
  barycoords = [1./2.*(1+sqrt(3./5.)) 1./2.*(1-sqrt(3./5.));
                1./2. 1./2.;
                1./2.*(1-sqrt(3./5.)) 1./2.*(1+sqrt(3./5.))]
  sum = 0.0
  for i = 1:3
    sum += weights[i] * f(barycoords[i,1]*p1 + barycoords[i,2]*p2)
  end
  return sum *length
end

"""
  setalldirichlet!(mesh, G, b)
Naively modifies in place the esm, G and the RHS, b for homogeneous Dirichlet BC.
"""
function setalldirichlet!(mesh::UniformTriangleMesh, G::Array{Float64, 2}, b::Array{Float64,2})
  # get boundary vertices
  V = reshape(1:size(mesh.vertices,1), mesh.m+1, mesh.n+1)
  V = flipdim(V,1) # make the matrix nodes look like picture
  interiorvertices = V[2:end-1, 2:end-1][:] # vector
  exteriorvertices = setdiff(1:size(mesh.vertices,1), interiorvertices)

  ind = trues(size(G)) # indices to zero out
  ind[interiorvertices,:] = false # keep interiorvertices rows
  for i in exteriorvertices # keep what's on the diagonal
    ind[i,i] = false
  end
  G[ind] = 0.0
  #= this is to set the Dirichlet boundary to some function.
  X,Y = ndgrid(linspace(0,1,mesh.n+1),linspace(0,1,mesh.m+1))
  g = X + Y
  bvals = g[exteriorvertices];
  =#
  b[exteriorvertices] = 0.0
end

end #MODULE
