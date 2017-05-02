%Finite Element CGP1 scheme on triangles 
%homework 2, problem 1

close all
clear all
tic
%define domain
x1 = 0;
x2 = 1;
y1 = 0;
y2 = 1;
n = 2^4;
%set up mesh
nx = n;
ny = n;
TriMesh = UnifTriMeshRectDomain(x1,x2,y1,y2,nx,ny);

%RHS function
RHS_funct = @(x,y) 2*pi*pi*sin(pi*x).*sin(pi*y);

%boundary conditions
%dirichlet on left and right boundaries of rectangular domain
dirichletconst = 0;
%neumann bcs on top and bottom boundaries of rectangular domain
Neumfunct_top =@(x,y) -pi*sin(pi*x).*cos(pi*y);
Neumfunct_bottom =@(x,y) pi*sin(pi*x).*cos(pi*y);

%array of all element stiffness matricies 
 AllESM = AllElmStiffnessArray_CGP1( TriMesh );
 %Assemble the global stiffness matrix, no boundary conditions are included
 %at this stage. 
 GSM = AssembleGSM_CGP1( TriMesh, AllESM );
 
 %array of local RHS used to assemble Global right hand side
 localRHS  = localRHS_CGP1( RHS_funct, TriMesh );
 %Assemble the global RHS
 GRHS = GlobalRHS_CGP1( localRHS, TriMesh );

 %now add neumann bcs to the RHS
 %GRHS_NBC = AddNeumannBCS_CGP1( TriMesh, GRHS, Neumfunct_top, Neumfunct_bottom );
 
 %add boundary conditions, first add dirichlet BCS
 [ GSM_DBC, GRHS_DBC ] = AddConstDirichlet_BCs_CGP1( GRHS, GSM, TriMesh, dirichletconst );

 
 %Solve the Linear System GSM_DBC(c)=GRHS_DBC, where c is a vector of
%constants such that sln = sum(ci*GlobalBasis_i)
sln = (GSM_DBC)\(GRHS_DBC);
slnMtx = reshape(sln, TriMesh.ny + 1, TriMesh.nx +1);

sln_max = max(sln)
toc
%Compute Exact Solution
%ExactSln = sin(pi.*RectMesh.x).*

ExactSln_Funct = @(x,y) sin(pi*x).*sin(pi*y);
ExactSln_partialx = @(x,y) pi*cos(pi*x).*sin(pi*y);
ExactSln_partialy = @(x,y) pi*sin(pi*x).*cos(pi*y);

ExactSln_grad = struct('partialx', ExactSln_partialx, 'partialy', ExactSln_partialy);

%evaluate the exact sln on the nodes of the mesh
X =TriMesh.meshgridXY(:,:,1);
Y =TriMesh.meshgridXY(:,:,2);
Z = ExactSln_Funct(X,Y);

exact_max = max(Z(:))

figure(3)
trisurf(TriMesh.elem_node, X,Y,Z)
title('The Exact Solution to the PDE')

figure(2)
trisurf(TriMesh.elem_node,X,Y,slnMtx)
title('Finite Element Solution to PDE')
xlabel('x')
ylabel('y')
rotate3d on;

%error_vec = Z(:)-sln;
%format long
%error = norm(error_vec, inf)

 
[ L2norm, H1norm ] = L2andH1errnorm( TriMesh, sln, ExactSln_Funct, ExactSln_grad)



