function [ localNeumannIntegrals ] = AddNeumannBCS_CGP1( TriMesh, GRHS, NeumFunct  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[ ~,top, bottom ] = Dirichlet_Neumann_BCs_rectdomain_CGP1_hw2( TriMesh );


localNeumannIntegrals = zeros(1,2,2*TriMesh.nx);
for i = 1: TriMesh.nx
    x1 = TriMesh.Xcoord(top(i));
    y1 = TriMesh.Ycoord(top(i));
    %confusion here. bc i need all three points of the element to find the
    %local basis functions for the integration along the neumann boundary
    %edge. So need to identify boundary elements, as well as the 
    phi1 = localbasis
integrand1 = 

end

