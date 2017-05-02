__precompile__()
module PoissonSolveTriP2
@everywhere push!(LOAD_PATH,"./")
include("UniformMesh.jl")

using UniformMesh
using Cubature # for quadrature
using PyPlot

export solve
export solveanddraw
export solveandnorm


"""
  poissonsolve(n)
Solves the problem and returns matrix of coeffs.
"""
function solve(n::Int64)
  # get mesh
  mesh = UniformTriangleMesh(n,n)
  numvertices = size(mesh.vertices,1)
  numedges = size(mesh.edges,1)
  # get global stiffness matrix
  G = gsm(mesh)
  # set RHS from quadrature over elements
  b = rhs(mesh, x -> 2*pi^2.*(sin(pi.*x[1]).*sin(pi.*x[2])))
  # set BC
  setalldirichlet!(mesh,G,b)
  # compute coeffs
  c = G\b
  return c
end


"""
  solveanddraw(n)
Solves the problem and plots the solution.
"""
function solveanddraw(n::Int64)
  c = solve(n)
  mesh = UniformTriangleMesh(n,n)
  A = [mesh.vertices c[1:(n+1)^2]]
  e1 = mesh.vertices[mesh.edges[:,1],:]
  e2 = mesh.vertices[mesh.edges[:,2],:]
  m = (mesh.vertices[mesh.edges[:,1],:] + mesh.vertices[mesh.edges[:,2],:])./2
  B = [m[:,1] m[:,2] c[(n+1)^2+1:end]]
  C = [A;B]
  #surf(A[:,1],A[:,2],c[1:(n+1)^2])
  #figure()
  #surf(C[:,1],C[:,2],C[:,3],cmap="viridis",edgecolor="None")

  fig = figure()
  ax = fig[:add_subplot](111, projection="3d")
  ax[:plot_trisurf](C[:,1],C[:,2],C[:,3], cmap="jet", edgecolors="black")
  ax[:set_title](string(n, "x", n), fontsize=22)
  ax[:set_xlabel]("X", fontsize=22)
  ax[:set_ylabel]("Y", fontsize=22)
  ax[:xaxis][:set_tick_params](labelsize=15)
  ax[:yaxis][:set_tick_params](labelsize=15)
  ax[:zaxis][:set_tick_params](labelsize=15)

end

"""
  solveandnorm(n)
Solves the problem for various grid sizes and computes the L2 and H1 norms for convergence purposes.
Prints a table.
"""
function solveandnorm()
  println(" ", "="^94)
  println("|", " "^4, "n", " "^4, "|", " "^4, "L2error", " "^5, "|", " "^5, "L2Conv",
          " "^5, "|", " "^5, "H1error", " "^4, "|", " "^5, "H1Conv", " "^5, "|",
          " "^3, "H1semierror", " "^2, "|")
  println(" ", "="^94)
  prevL2error = 1.0
  prevH1error = 1.0
  nold = 8
  for n in [2, 4, 8, 16, 32]
    c = solve(n)
    errL2 = errorL2(c,n)
    errH1 = errorH1(c,n)
    convergenceL2 = log(2, prevL2error/errL2)
    prevL2error = errL2
    convergenceH1 = log(2, prevH1error/sqrt(errL2^2 + errH1^2))
    prevH1error = sqrt(errL2^2 + errH1^2)
    print("|", " "^3)
    @printf("%3d", n)
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
    println(" "^3, "|")
  end
  println(" ", "="^94)
end

"""
  errorL2()
Compute the L2 error between the exact solution and the approximate solution.
TODO needs to be over 6 points
"""
function errorL2(c::Array{Float64, 2}, n::Int64)
  mesh = UniformTriangleMesh(n,n)
  error = 0.0
  for t = 1:size(mesh.triangles,1)
    trierror = 0.0
    p1 = mesh.vertices[mesh.triangles[t,1],:]
    p2 = mesh.vertices[mesh.triangles[t,2],:]
    p3 = mesh.vertices[mesh.triangles[t,3],:]
    numvertices = size(mesh.vertices,1)
    edgeindices = triangleedgeindices(mesh, mesh.triangles[t,:])
    edgeindices += numvertices # edge numbering starts after vertices
    u(x) = sin(pi*x[1])*sin(pi*x[2])
    function uh()
      c1,c2,c3 = c[mesh.triangles[t,:]]
      c4,c5,c6 = c[edgeindices]
      return [c1,c2,c3,c4,c5,c6]
    end
    # integrate the difference of the exact and approximate solutions
    trierror = trigaussquad(x -> norm(u(x) - uh()[1] * φ1(x,p1,p2,p3)
                                           - uh()[2] * φ2(x,p1,p2,p3)
                                           - uh()[3] * φ3(x,p1,p2,p3)
                                           - uh()[4] * φ23(x,p1,p2,p3)
                                           - uh()[5] * φ13(x,p1,p2,p3)
                                           - uh()[6] * φ12(x,p1,p2,p3) )^2, p1, p2, p3)

    error += trierror
  end
  return sqrt(error)
end

"""
  errorH1()
Compute the H1 error between the exact solution and the approximate solution.
TODO needs to be over 6 points.
"""
function errorH1(c::Array{Float64, 2}, n::Int64)
  mesh = UniformTriangleMesh(n,n)
  error = 0.0
  for t = 1:size(mesh.triangles,1)
    trierror = 0.0
    p1 = mesh.vertices[mesh.triangles[t,1],:]
    p2 = mesh.vertices[mesh.triangles[t,2],:]
    p3 = mesh.vertices[mesh.triangles[t,3],:]
    numvertices = size(mesh.vertices,1)
    edgeindices = triangleedgeindices(mesh, mesh.triangles[t,:])
    edgeindices += numvertices # edge numbering starts after vertices
    u(x) = [pi*cos(pi*x[1])*sin(pi*x[2]),pi*sin(pi*x[1])*cos(pi*x[2])]
    function uh()
      c1,c2,c3 = c[mesh.triangles[t,:]]
      c4,c5,c6 = c[edgeindices]
      return [c1,c2,c3,c4,c5,c6]
    end
    # integrate the difference of the exact and approximate solution
    trierror = trigaussquad(x -> norm(u(x) - uh()[1] .* ∇φ1(x,p1,p2,p3)
                                           - uh()[2] .* ∇φ2(x,p1,p2,p3)
                                           - uh()[3] .* ∇φ3(x,p1,p2,p3)
                                           - uh()[4] .* ∇φ23(x,p1,p2,p3)
                                           - uh()[5] .* ∇φ13(x,p1,p2,p3)
                                           - uh()[6] .* ∇φ12(x,p1,p2,p3))^2, p1, p2, p3)
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
Do quadrature over a triangular element for RHS us P2 basis functions.
"""
function elementrhs(p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1},
                    m1::Array{Float64,1}, m2::Array{Float64,1}, m3::Array{Float64,1}, f::Function)
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  b1 = trigaussquad(x -> (f(x) * φ1(x,p1,p2,p3)), p1, p2, p3)
  b2 = trigaussquad(x -> (f(x) * φ2(x,p1,p2,p3)), p1, p2, p3)
  b3 = trigaussquad(x -> (f(x) * φ3(x,p1,p2,p3)), p1, p2, p3)
  b4 = trigaussquad(x -> (f(x) * φ23(x,p1,p2,p3)), p1, p2, p3)
  b5 = trigaussquad(x -> (f(x) * φ13(x,p1,p2,p3)), p1, p2, p3)
  b6 = trigaussquad(x -> (f(x) * φ12(x,p1,p2,p3)), p1, p2, p3)
  return [b1; b2; b3; b4; b5; b6;]
end

"""
  rhs(mesh)
Assemble the vector b (RHS) by summing over elements.
"""
function rhs(mesh::UniformTriangleMesh, f::Function)
  b = zeros(Float64, size(mesh.vertices,1) + size(mesh.edges,1), 1) # preallocate

  for i = 1:size(mesh.triangles,1)
    numvertices = size(mesh.vertices,1)
    numedges = size(mesh.edges,1)
    midpoints = edgestomidpoints(mesh, triangletoedges(mesh.triangles[i,:]))
    elemb = elementrhs(mesh.vertices[mesh.triangles[i,1],:],
                       mesh.vertices[mesh.triangles[i,2],:],
                       mesh.vertices[mesh.triangles[i,3],:],
                       midpoints[1,:],
                       midpoints[2,:],
                       midpoints[3,:],
                       f)
    triangleindices = mesh.triangles[i,:]
    edgeindices = triangleedgeindices(mesh, mesh.triangles[i,:])
    edgeindices += numvertices # edge numbering starts after vertices
    b[triangleindices] += elemb[1:3]
    b[edgeindices] += elemb[4:6]
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
  E11 = trigaussquad(x -> dot(∇φ1(x,p1,p2,p3) , ∇φ1(x,p1,p2,p3)), p1, p2, p3)
  E12 = trigaussquad(x -> dot(∇φ1(x,p1,p2,p3) , ∇φ2(x,p1,p2,p3)), p1, p2, p3)
  E13 = trigaussquad(x -> dot(∇φ1(x,p1,p2,p3) , ∇φ3(x,p1,p2,p3)), p1, p2, p3)
  E14 = trigaussquad(x -> dot(∇φ1(x,p1,p2,p3) , ∇φ23(x,p1,p2,p3)), p1, p2, p3)
  E15 = trigaussquad(x -> dot(∇φ1(x,p1,p2,p3) , ∇φ13(x,p1,p2,p3)), p1, p2, p3)
  E16 = trigaussquad(x -> dot(∇φ1(x,p1,p2,p3) , ∇φ12(x,p1,p2,p3)), p1, p2, p3)

  E22 = trigaussquad(x -> dot(∇φ2(x,p1,p2,p3) , ∇φ2(x,p1,p2,p3)), p1, p2, p3)
  E23 = trigaussquad(x -> dot(∇φ2(x,p1,p2,p3) , ∇φ3(x,p1,p2,p3)), p1, p2, p3)
  E24 = trigaussquad(x -> dot(∇φ2(x,p1,p2,p3) , ∇φ23(x,p1,p2,p3)), p1, p2, p3)
  E25 = trigaussquad(x -> dot(∇φ2(x,p1,p2,p3) , ∇φ13(x,p1,p2,p3)), p1, p2, p3)
  E26 = trigaussquad(x -> dot(∇φ2(x,p1,p2,p3) , ∇φ12(x,p1,p2,p3)), p1, p2, p3)

  E33 = trigaussquad(x -> dot(∇φ3(x,p1,p2,p3) , ∇φ3(x,p1,p2,p3)), p1, p2, p3)
  E34 = trigaussquad(x -> dot(∇φ3(x,p1,p2,p3) , ∇φ23(x,p1,p2,p3)), p1, p2, p3)
  E35 = trigaussquad(x -> dot(∇φ3(x,p1,p2,p3) , ∇φ13(x,p1,p2,p3)), p1, p2, p3)
  E36 = trigaussquad(x -> dot(∇φ3(x,p1,p2,p3) , ∇φ12(x,p1,p2,p3)), p1, p2, p3)

  E44 = trigaussquad(x -> dot(∇φ23(x,p1,p2,p3) , ∇φ23(x,p1,p2,p3)), p1, p2, p3)
  E45 = trigaussquad(x -> dot(∇φ23(x,p1,p2,p3) , ∇φ13(x,p1,p2,p3)), p1, p2, p3)
  E46 = trigaussquad(x -> dot(∇φ23(x,p1,p2,p3) , ∇φ12(x,p1,p2,p3)), p1, p2, p3)

  E55 = trigaussquad(x -> dot(∇φ13(x,p1,p2,p3) , ∇φ13(x,p1,p2,p3)), p1, p2, p3)
  E56 = trigaussquad(x -> dot(∇φ13(x,p1,p2,p3) , ∇φ12(x,p1,p2,p3)), p1, p2, p3)

  E66 = trigaussquad(x -> dot(∇φ12(x,p1,p2,p3) , ∇φ12(x,p1,p2,p3)), p1, p2, p3)

  E = zeros(Float64, 6, 6)

      #  vertices    edges
  E = [E11 E12 E13 E14 E15 E16;
       E12 E22 E23 E24 E25 E26; # vertices
       E13 E23 E33 E34 E35 E36;
       E14 E24 E34 E44 E45 E46;
       E15 E25 E35 E45 E55 E56; # edges
       E16 E26 E36 E46 E56 E66]
end

"""
  gsm(mesh)
Assemble the gsm from all of the element stiffness matrices.
"""
function gsm(mesh::UniformTriangleMesh)
  numvertices = size(mesh.vertices,1)
  numedges = size(mesh.edges,1)
  G = zeros(Float64, numvertices + numedges, numvertices + numedges) # preallocate
  for i = 1:size(mesh.triangles,1)
    E = esmtri(mesh.vertices[mesh.triangles[i,1],:],
               mesh.vertices[mesh.triangles[i,2],:],
               mesh.vertices[mesh.triangles[i,3],:])
    triangleindices = mesh.triangles[i,:]
    edgeindices = triangleedgeindices(mesh, mesh.triangles[i,:])
    edgeindices += numvertices # edge numbering starts after vertices
    G[triangleindices, triangleindices] += E[1:3,1:3]
    G[triangleindices, edgeindices] += E[1:3,4:6]
    G[edgeindices, triangleindices] += E[4:6,1:3]
    G[edgeindices, edgeindices] += E[4:6,4:6]
  end

  return G
end


"""
  setalldirichlet!(mesh, G, b)
Naively modifies in place the esm, G and the RHS, b for homogeneous Dirichlet BC.
"""
function setalldirichlet!(mesh::UniformTriangleMesh, G::Array{Float64, 2}, b::Array{Float64,2})
  numvertices = size(mesh.vertices,1)
  # get boundary vertices
  V = reshape(1:numvertices, mesh.m+1, mesh.n+1)
  V = flipdim(V,1) # make the matrix nodes look like picture
  interiorvertices = V[2:end-1, 2:end-1][:] # vector
  exteriorvertices = setdiff(1:numvertices, interiorvertices)
  ind = trues((mesh.m + 1)^2, (mesh.n + 1)^2) # indices to zero out
  ind[interiorvertices,:] = false # keep interiorvertices rows
  for i in exteriorvertices # keep what's on the diagonal
    ind[i,i] = false
  end
  B = G[1:(mesh.m + 1)^2,1:(mesh.n + 1)^2]
  B[ind] = 0.0
  G[1:(mesh.m + 1)^2,1:(mesh.n + 1)^2] = B
  B = G[exteriorvertices,numvertices+1:end] = 0.0
  b[exteriorvertices] = 0.0 # g(exteriorvertices)
  setedgesdirichlet!(mesh,G,b,exteriorvertices)
end

function setedgesdirichlet!(mesh::UniformTriangleMesh, G::Array{Float64, 2}, b::Array{Float64,2},exteriorvertices::Array{Int64,1})
  # if both vertices of an edge are in the boundary, then set corresponding row in G,b to zero
  numvertices = size(mesh.vertices,1)
  numedges = size(mesh.edges,1)
  for i = 1:numedges
    index1 = find(exteriorvertices .== mesh.edges[i,:][1])
    index2 = find(exteriorvertices .== mesh.edges[i,:][2])
    if isempty(index1) || isempty(index2)
      continue
    else
      G[i+numvertices,1:i+numvertices-1] = 0.0
      G[i+numvertices,i+numvertices+1:end] = 0.0
      b[i+numvertices] = 0.0
    end
  end
end


# Barycentric coords
function λ1(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  0.5 * abs(det([1 x[1] x[2]; 1 p2[1] p2[2]; 1 p3[1] p3[2]]))/area
end

function λ2(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  0.5 * abs(det([1 p1[1] p1[2]; 1 x[1] x[2]; 1 p3[1] p3[2]]))/area
end

function λ3(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  0.5 * abs(det([1 p1[1] p1[2]; 1 p2[1] p2[2]; 1 x[1] x[2]]))/area
end

# Barycentric coords gradients
function ∇λ1(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  [(p2[2] - p3[2]),(p3[1] - p2[1])]/(2*area)
end

function ∇λ2(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  [(p3[2] - p1[2]),(p1[1]) - p3[1]]/(2*area)
end

function ∇λ3(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  [(p1[2] - p2[2]),(p2[1] - p1[1])]/(2*area)
end

# P2 Basis functions
φ1(x::Array{Float64,1} ,p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1}) =
  λ1(x,p1,p2,p3)*(2*λ1(x,p1,p2,p3) - 1)

φ2(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1}) =
  λ2(x,p1,p2,p3)*(2*λ2(x,p1,p2,p3) - 1)

φ3(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1}) =
  λ3(x,p1,p2,p3)*(2*λ3(x,p1,p2,p3) - 1)

φ12(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1}) =
  4*λ1(x,p1,p2,p3)*λ2(x,p1,p2,p3)

φ13(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1}) =
  4*λ1(x,p1,p2,p3)*λ3(x,p1,p2,p3)

φ23(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1}) =
  4*λ2(x,p1,p2,p3)*λ3(x,p1,p2,p3)

# P2 Basis function gradients
function ∇φ1(x::Array{Float64,1} ,p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * det([1 p1[1] p1[2]; 1 p2[1] p2[2]; 1 p3[1] p3[2]])
  ∇λ1(x,p1,p2,p3)*(4*λ1(x,p1,p2,p3) - 1)
end

function ∇φ2(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * det([1 p1[1] p1[2]; 1 p2[1] p2[2]; 1 p3[1] p3[2]])
  ∇λ2(x,p1,p2,p3)*(4*λ2(x,p1,p2,p3) - 1)
end

function ∇φ3(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * det([1 p1[1] p1[2]; 1 p2[1] p2[2]; 1 p3[1] p3[2]])
  ∇λ3(x,p1,p2,p3)*(4*λ3(x,p1,p2,p3) - 1)
end

function ∇φ12(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * det([1 p1[1] p1[2]; 1 p2[1] p2[2]; 1 p3[1] p3[2]])
  4*λ1(x,p1,p2,p3)*∇λ2(x,p1,p2,p3) + 4*λ2(x,p1,p2,p3)*∇λ1(x,p1,p2,p3)
end

function ∇φ13(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * det([1 p1[1] p1[2]; 1 p2[1] p2[2]; 1 p3[1] p3[2]])
  4*λ1(x,p1,p2,p3)*∇λ3(x,p1,p2,p3) + 4*λ3(x,p1,p2,p3)*∇λ1(x,p1,p2,p3)
end

function ∇φ23(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * det([1 p1[1] p1[2]; 1 p2[1] p2[2]; 1 p3[1] p3[2]])
  4*λ2(x,p1,p2,p3)*∇λ3(x,p1,p2,p3) + 4*λ3(x,p1,p2,p3)*∇λ2(x,p1,p2,p3)
end
#=
function ∇φ1(x::Array{Float64,1} ,p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  [-(p2[2] - p3[2])*(-4*p2[1]*p3[2] + 4*p2[1]*x[2] + 4*p2[2]*(p3[1] - x[1]) + 4*p3[2]*x[1] - 4*p3[1]*x[2] + 1),
   -(p3[1] - p2[1])*(-4*p2[1]*p3[2] + 4*p2[1]*x[2] + 4*p2[2]*(p3[1] - x[1]) + 4*p3[2]*x[1] - 4*p3[1]*x[2] + 1)]./(2*area)
end

function ∇φ2(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  [-(p1[2] - p3[2])*(-4*p1[1]*p3[2] + 4*p1[1]*x[2] + 4*p1[2]*(p3[1] - x[1]) + 4*p3[2]*x[1] - 4*p3[1]*x[2] - 1),
   -(p1[1] - p3[1])*(4*p1[1]*p3[2] - 4*p1[1]*x[2] + 4*p1[2]*x[1] - 4*p3[1]*p1[2] - 4*p3[2]*x[1] + 4*p3[1]*x[2] + 1)]./(2*area)
 end

function ∇φ3(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  [-(p1[2] - p2[2])*(-4*p1[1]*p2[2] + 4*p1[1]*x[2] + 4*p1[2]*(p2[1] - x[1]) - 4*p2[1]*x[2] + 4*p2[2]*x[1] + 1),
  -(p1[1] - p2[1])*(4*p1[1]*p2[2] - 4*p1[1]*x[2] - 4*p1[2]*p2[1] + 4*p1[2]*x[1] + 4*p2[1]*x[2] - 4*p2[2]*x[1] - 1)]./(2*area)
end

function ∇φ12(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  [4*(p2[2] - p3[2])*(p1[1]*(x[2] - p3[2]) + p1[2]*(p3[1] - x[1]) - p3[1]*x[2] + p3[2]*x[1]) + 4*(p1[2] - p3[2])*(p2[1]*(x[2] - p3[2]) + p2[2]*(p3[1] - x[1]) - p3[1]*x[2] + p3[2]*x[1]),
   4*(p3[1] - p2[1])*(p1[1]*(x[2] - p3[2]) + p1[2]*(p3[1] - x[1]) - p3[1]*x[2] + p3[2]*x[1]) + 4*(p1[1] - p3[1])*(p2[1]*(p3[2] - x[2]) + p2[2]*(x[1] - p3[1]) + p3[1]*x[2] - p3[2]*x[1])]./(2*area)
 end

function ∇φ13(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  [4*(p2[2] - p3[2])*(p1[1]*(p2[2] - x[2]) + p1[2]*(x[1] - p2[1]) + p2[1]*x[2] - p2[2]*x[1]) + 4*(p1[2] - p2[2])*(p2[1]*(p3[2] - x[2]) + p2[2]*(x[1] - p3[1]) + p3[1]*x[2] - p3[2]*x[1]),
   4*(p2[1] - p3[1])*(p1[1]*(x[2] - p2[2]) + p1[2]*(p2[1] - x[1]) - p2[1]*x[2] + p2[2]*x[1]) + 4*(p2[1] - p1[1])*(p2[1]*(p3[2] - x[2]) + p2[2]*(x[1] - p3[1]) + p3[1]*x[2] - p3[2]*x[1])]./(2*area)
 end

function ∇φ23(x::Array{Float64,1}, p1::Array{Float64,1}, p2::Array{Float64,1}, p3::Array{Float64,1})
  area = 0.5 * abs(det([p1[1] p1[2] 1; p2[1] p2[2] 1; p3[1] p3[2] 1]))
  [4*(p1[2] - p3[2])*(p1[1]*(x[2] - p2[2]) + p1[2]*(p2[1] - x[1]) - p2[1]*x[2] + p2[2]*x[1]) + 4*(p1[2] - p2[2])*(p1[1]*(x[2] - p3[2]) + p1[2]*(p3[1] - x[1]) - p3[1]*x[2] + p3[2]*x[1]),
   4*(p1[1] - p3[1])*(p1[1]*(p2[2] - x[2]) + p1[2]*(x[2] - p2[1]) + p2[1]*x[2] - p2[2]*x[1]) + 4*(p1[1] - p2[1])*(p1[1]*(p3[2] - x[2]) + p1[2]*(x[1] - p3[1]) + p3[1]*x[2] - p3[2]*x[1])]./(2*area)
 end
=#

function ndgrid{T}(v1::AbstractArray{T}, v2::AbstractArray{T})
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repmat(v1, 1, n), repmat(v2, m, 1))
end

end # module
