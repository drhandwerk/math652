function [ localRHS ] = localRHS_CGP1( RHS_funct, TriMesh )
%UNTITLED3 Summary of this function goes here
%output: an array consisting of local rhs vectors

localRHS = zeros(1,3,TriMesh.numbelem);

for i = 1:TriMesh.numbelem
    NI = TriMesh.elem_node(i,:);
        x1 = TriMesh.Xcoord(NI(1));
        y1 = TriMesh.Ycoord(NI(1));
        x2 = TriMesh.Xcoord(NI(2));
        y2 = TriMesh.Ycoord(NI(2));
        x3 = TriMesh.Xcoord(NI(3));
        y3 = TriMesh.Ycoord(NI(3));
    localbasis = localbasisfunct_CGP1( x1,y1,x2,y2,x3,y3 );
    
    integrand1 = @(x,y)localbasis.phi1(x,y).*RHS_funct(x,y);
    integrand2 = @(x,y)localbasis.phi2(x,y).*RHS_funct(x,y);
    integrand3 = @(x,y)localbasis.phi3(x,y).*RHS_funct(x,y);
    
    localRHS(1,1,i) = GaussQuad_Tri3pt( x1,y1, x2,y2,x3,y3, integrand1 );
    localRHS(1,2,i) = GaussQuad_Tri3pt( x1,y1, x2,y2,x3,y3, integrand2 );
    localRHS(1,3,i) = GaussQuad_Tri3pt( x1,y1, x2,y2,x3,y3, integrand3 );
    
end

end

