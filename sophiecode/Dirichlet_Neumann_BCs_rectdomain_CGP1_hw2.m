function [ leftandright,top, bottom ] = Dirichlet_Neumann_BCs_rectdomain_CGP1_hw2( TriMesh )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

leftbndry = TriMesh.bndrynodes(1:TriMesh.ny +1);

rightbndry = TriMesh.bndrynodes(TriMesh.ny + 2*TriMesh.nx : end);

leftandright = vertcat(leftbndry, rightbndry);

%topandbottom = setdiff(TriMesh.bndrynodes, leftandright);
%top = topandbottom(1:2:end);
%bottom = topandbottom(2:2:end);

top = TriMesh.nodeindex(1,:);
bottom = TriMesh.nodeindex(end,:);

end

