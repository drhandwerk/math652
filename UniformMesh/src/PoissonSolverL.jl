__precompile__()
module PoissonSolverL
@everywhere push!(LOAD_PATH,"./")
include("UniformMesh.jl")

using UniformMesh
using Cubature # for quadrature
using PyPlot

export poissonsolve
export solveanddraw
export solveandnorm

"""
 atan2custom(y,x)
So no discontinuity at pi
"""
function atan2custom(y::Array{Float64,1},x::Array{Float64,1})
  angles = atan2(y,x)
  angles[angles.<0] += 2pi
  return angles
end
function atan2custom(y::Float64,x::Float64)
  atan2(y,x) > 0 ? atan2(y,x) : atan2(y,x) + 2pi
end

"""
  poissonsolve(n)
Solves the problem and returns matrix of coeffs.
"""
function poissonsolve(n::Float64)
  # get mesh
  mesh = DistMeshTriangleMesh(n)
  # get global stiffness matrix
  G = gsm(mesh)
  # set RHS from quadrature over elements
  b = rhs(mesh, x -> 0.0)
  # set BC
  #setneumann!(mesh,b)

  setalldirichlet!(mesh,G,b,(x,y) -> (x.^2+y.^2).^(1/3).*sin(2/3.*atan2custom(y,x)))
  # compute coeffs
  c = G\b
  # return as matrix
  return c
end



"""
  solveanddraw(n)
Solves the problem and plots the solution.
"""
function solveanddraw(n::Float64)
  c = poissonsolve(n)
  trimesh = DistMeshTriangleMesh(n)
  fig = figure()
  ax = fig[:add_subplot](111, projection="3d")
  ax[:plot_trisurf](trimesh.vertices[:,1],trimesh.vertices[:,2], triangles=trimesh.triangles-1, c[:], alpha=1, cmap="viridis", edgecolors="None")
  ax[:set_title](string("L-Shape"), fontsize=18)
  ax[:set_xlabel]("X", fontsize=16)
  ax[:set_ylabel]("Y", fontsize=16)
  ax[:xaxis][:set_tick_params](labelsize=8)
  ax[:yaxis][:set_tick_params](labelsize=8)
  ax[:zaxis][:set_tick_params](labelsize=8)
end

"""
  solveandtime()
Solve and time the problem for various mesh sizes. The mesh sizes specified must have already been
created using DistMesh in MATLAB.
"""
function solveandtime()
  println(" ", "="^30)
  println("|", " "^5, "h0", " "^6, "|", " "^7, "Time", " "^5, "|")
  println(" ", "="^30)

  for h0 in [.4, .2, .1, .05, .025]
    time = @elapsed c = poissonsolve(h0)
    print("|", " "^3)
    @printf("%2.4f", h0)
    print(" "^3, "|", " "^3)
    @printf("%6.8f", time)
    println(" "^3, "|")
  end
  println(" ", "="^30)
end
"""
  solveandnorm()
Solves the problem for various grid sizes and computes the L2 and H1 norms for convergence purposes.
Prints a table.
"""
function solveandnorm()
  println(" ", "="^114)
  println("|", " "^4, "h0", " "^5, "|", " "^4, "L2error", " "^5, "|", " "^5, "L2Conv",
          " "^5, "|", " "^5, "H1error", " "^4, "|", " "^5, "H1Conv", " "^5, "|",
          " "^3, "H1semierror", " "^2, "|", " "^7, "Time", " "^5, "|")
  println(" ", "="^114)
  prevL2error = 1.0
  prevH1error = 1.0
  for h0 in [.4, .2, .1, .05, .025]
    time = @elapsed c = poissonsolve(h0)
    function atansign(y::Float64,x::Float64)
      atan2(y,x) < 0 ? -1 : 1
    end
    errL2 = errorL2(c, h0, x -> (x[1].^2+x[2].^2).^(1/3).*sin(2./3.*atan2custom(x[2],x[1])))
    errH1 = errorH1(c, h0, x -> [(cos(atan2custom(x[2],x[1])).*2./3.*sqrt(x[2].^2+x[1].^2).^(-1./3.).*sin(2./3.*atan2custom(x[2],x[1]))
                                 - 2./3.*sin(atan2custom(x[2],x[1])).*sqrt(x[2].^2+x[1].^2).^(-1./3.).*cos(2./3.*atan2custom(x[2],x[1]))),
                                 (sin(atan2custom(x[2],x[1])).*2./3.*sqrt(x[2].^2+x[1].^2).^(-1./3.).*sin(2./3.*atan2custom(x[2],x[1]))
                                 + 2./3.*cos(atan2custom(x[2],x[1])).*sqrt(x[2].^2+x[1].^2).^(-1./3.).*cos(2./3.*atan2custom(x[2],x[1])))])
    convergenceL2 = log(2, prevL2error/errL2)
    prevL2error = errL2
    convergenceH1 = log(2, prevH1error/sqrt(errL2^2 + errH1^2))
    prevH1error = sqrt(errL2^2 + errH1^2)
    print("|", " "^3)
    @printf("%2.3f", h0)
    print(" "^3, "|", " "^3)
    @printf("%6.8f", errL2)
    print(" "^3, "|", " "^3)
    @printf("%6.8f", convergenceL2)
    print(" "^3, "|", " "^3)
    @printf("%6.8f", sqrt(errL2^2 + errH1^2))
    print(" "^3, "|", " "^3)
    @printf("%6.8f", convergenceH1)
    print(" "^3, "|", " "^3)
    @printf("%6.8f", errH1)
    print(" "^3, "|", " "^3)
    @printf("%6.8f", time)
    println(" "^3, "|")
  end
  println(" ", "="^114)
end

"""
  errorL2(c, h0, u)
Compute the L2 error between the exact solution 'u' and the approximate solution 'c'.
The initial edge length 'h0' is from DistMesh in MATLAB.
"""
function errorL2(c::Array{Float64, 2}, h0::Float64, u::Function)
  trimesh = DistMeshTriangleMesh(h0)
  error = 0.0
  for t = 1:size(trimesh.triangles,1)
    trierror = 0.0
    p1 = trimesh.vertices[trimesh.triangles[t,1],:]
    p2 = trimesh.vertices[trimesh.triangles[t,2],:]
    p3 = trimesh.vertices[trimesh.triangles[t,3],:]
    area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
    c1,c2,c3 = c[:][trimesh.triangles[t,:]]
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
  errorH1()
Compute the H1 error between the exact solution and the approximate solution.
"""
function errorH1(c::Array{Float64, 2}, h0::Float64, Du::Function)
  trimesh = DistMeshTriangleMesh(h0)
  error = 0.0
  for t = 1:size(trimesh.triangles,1)
    trierror = 0.0
    p1 = trimesh.vertices[trimesh.triangles[t,1],:]
    p2 = trimesh.vertices[trimesh.triangles[t,2],:]
    p3 = trimesh.vertices[trimesh.triangles[t,3],:]
    area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
    c1,c2,c3 = c[:][trimesh.triangles[t,:]] # approximate solution
    uh(x) = x == p1 ? c1 : x == p2 ? c2 : c3
    # integrate the difference of the exact and approximate solution
    trierror = trigaussquad(x -> norm(Du(x) - uh(p1) .* [p2[2]-p3[2],p3[1]-p2[1]]./(2.0*area)
                                            - uh(p2) .* [p3[2]-p1[2],p1[1]-p3[1]]./(2.0*area)
                                            - uh(p3) .* [p1[2]-p2[2],p2[1]-p1[1]]./(2.0*area))^2, p1, p2, p3)
    error += trierror
  end
  return sqrt(error)
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
  elementrhs()
Do quadrature over a triangular element for RHS.
"""
function elementrhs(p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1}, f::Function)
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  b1 = trigaussquad(x -> f(x) * 0.5 * abs(det([x[1] x[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))/area, p1, p2, p3)
  b2 = trigaussquad(x -> f(x) * 0.5 * abs(det([p1[1] p1[2] 1; x[1] x[2] 1; p3[1] p3[2] 1]))/area, p1, p2, p3)
  b3 = trigaussquad(x -> f(x) * 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; x[1] x[2] 1]))/area, p1, p2, p3)
  return [b1; b2; b3]
end

"""
  rhs(mesh)
Assemble the vector b (RHS) by summing over elements.
"""
function rhs(mesh::DistMeshTriangleMesh, f::Function)
  b = zeros(Float64, size(mesh.vertices,1), 1) # preallocate
  for i = 1:size(mesh.triangles,1)
    elemb = elementrhs(mesh.vertices[mesh.triangles[i,1],:],
                       mesh.vertices[mesh.triangles[i,2],:],
                       mesh.vertices[mesh.triangles[i,3],:],
                       f)
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
function gsm(mesh::DistMeshTriangleMesh)
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
Modify b using quadrature on the boundary. Don't actually need to send G in here.
"""
function setneumann!(mesh::DistMeshTriangleMesh, b::Array{Float64, 2})
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
    #p1 = mesh.vertices[mesh.triangles[e,1],:]
    p2 = mesh.vertices[mesh.triangles[e,2],:]
    p3 = mesh.vertices[mesh.triangles[e,3],:]
    length = sqrt((p2[1]-p3[1])^2 + (p2[2]-p3[2])^2)
    b1 = lingaussquad(x -> pi*sin(pi*x[1])*cos(pi*x[2]) * sqrt((p2[1]-x[1])^2 + (p2[2]-x[2])^2)/length, p3, p2)
    b2 = lingaussquad(x -> pi*sin(pi*x[1])*cos(pi*x[2]) * sqrt((p3[1]-x[1])^2 + (p3[2]-x[2])^2)/length, p3, p2)
    b[mesh.triangles[e,:]] += [0.0; b2; b1] #need to add instead of subtract because of how nodes/G are ordered/positioned

  end
  for e in bottomelements
    p1 = mesh.vertices[mesh.triangles[e,1],:]
    p2 = mesh.vertices[mesh.triangles[e,2],:]
    #p3 = mesh.vertices[mesh.triangles[e,3],:]
    length = sqrt((p2[1]-p1[1])^2 + (p2[2]-p1[2])^2)
    b1 = lingaussquad(x -> -pi*sin(pi*x[1])*cos(pi*x[2]) * sqrt((p2[1]-x[1])^2 + (p2[2]-x[2])^2)/length, p1, p2)
    b2 = lingaussquad(x -> -pi*sin(pi*x[1])*cos(pi*x[2]) * sqrt((p1[1]-x[1])^2 + (p1[2]-x[2])^2)/length, p1, p2)
    b[mesh.triangles[e,:]] += [b1; b2; 0.0]

  end
end

"""
  setdirichlet!(mesh, G, b)
Naively modifies in place the esm, G and the RHS, b for homogeneous Dirichlet BC on left and right.
"""
function setdirichlet!(mesh::DistMeshTriangleMesh, G::Array{Float64, 2}, b::Array{Float64,2})
  # get all boundary vertices
  V = reshape(1:size(mesh.vertices,1), mesh.m+1, mesh.n+1)
  V = flipdim(V,1)
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
    ind[i,i] = false
  end
  for i in rightvertices
    ind[i,i] = false
  end
  G[ind] = 0.0
  #= this is to set the Dirichlet boundary to some function.
  X,Y = ndgrid(linspace(0,1,mesh.n+1),linspace(0,1,mesh.m+1))
  g = X + Y
  bvals = g[exteriorvertices];
  =#
  b[leftvertices] = 0.0
  b[rightvertices] = 0.0
end

"""
  setalldirichlet!(mesh, G, b, g)
Naively modifies in place the esm, G and the RHS, b for homogeneous Dirichlet BC.
Set the boundary nodes to the function 'g(x,y)'.
"""
function setalldirichlet!(mesh::DistMeshTriangleMesh, G::Array{Float64, 2}, b::Array{Float64,2}, g::Function)
  ind = falses(size(G)) # indices to zero out
  for i in mesh.boundaryvertices # keep what's on the diagonal
    ind[i,1:i-1] = true
    ind[i,i+1:end] = true
    G[i,i] = 1.0 # for zero case doesn't matter. for nonzero case don't want to rescale rhs
  end
  G[ind] = 0.0
  bvals = g(mesh.vertices[:,1],mesh.vertices[:,2]);
  b[mesh.boundaryvertices] = bvals[mesh.boundaryvertices]
end

function ndgrid{T}(v1::AbstractArray{T}, v2::AbstractArray{T})
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repmat(v1, 1, n), repmat(v2, m, 1))
end

end # module
