function [ GRHS] = AddNeumannBCS_CGP1( TriMesh, GRHS, Neumfunct_top, Neumfunct_bottom )
%computing local neumann boundary conditions and assembling the GRHS to
%incorporate these neumann bcs. 

%top boundary for Neumann

for i =1:size(TriMesh.topbndryelem,1)
    NI = TriMesh.topbndryelem(i,:);
    x1 = TriMesh.Xcoord(NI(1));
    y1 = TriMesh.Ycoord(NI(1));
    x2 = TriMesh.Xcoord(NI(2));
    y2 = TriMesh.Ycoord(NI(2));
    x3 = TriMesh.Xcoord(NI(3));
    y3 = TriMesh.Ycoord(NI(3));
    
   localbasisfunct = localbasisfunct_CGP1(x1,y1,x2,y2,x3,y3);
    
   %phi1_other = @(x,y) sqrt((x2-x).^2+(y2-y).^2);
   %phi2_other = @(x,y) sqrt((x-x1).^2+(y-y1).^2);
   
    integrand1 = @(x,y)localbasisfunct.phi1(x,y).*Neumfunct_top(x,y);
    integrand2 = @(x,y)localbasisfunct.phi2(x,y).*Neumfunct_top(x,y);
    
    %integrand1 = @(x,y)phi1_other(x,y).*Neumfunct_top(x,y);
    %integrand2 = @(x,y)phi2_other(x,y).*Neumfunct_top(x,y);
    
    integrand1x = @(x)integrand1(x,1);
    integrand2x = @(x)integrand2(x,1);
    
    NeumIntegral1 = GaussQuad_Line3pt(x2,x1,integrand1x);
    NeumIntegral2 = GaussQuad_Line3pt(x2,x1,integrand2x);
    
    GRHS(NI(1)) = GRHS(NI(1))-NeumIntegral1;
    GRHS(NI(2)) = GRHS(NI(2))-NeumIntegral2;
end

%the following is bc the corners are associated with the dirichlet boundary
%condition, so to avoid double assigning boundary conditions at that node,
%we do the the first and last elements separately. 
%first element-- the node index (2) is not included in neumann
%THIS MAY NOT BE NECESSARY
%     NIone = TriMesh.topbndryelem(1,:);
%     x1 = TriMesh.Xcoord(NIone(1));
%     y1 = TriMesh.Ycoord(NIone(1));
%     x2 = TriMesh.Xcoord(NIone(2));
%     y2 = TriMesh.Ycoord(NIone(2));
%     x3 = TriMesh.Xcoord(NIone(3));
%     y3 = TriMesh.Ycoord(NIone(3));
%     
%    localbasisfunct = localbasisfunct_CGP1(x1,y1,x2,y2,x3,y3);
%     
%     integrand1 = @(x,y,n)localbasisfunct.phi1(x,y).*Neumfunct(x,y,n);
%     
%     integrand1x = @(x)integrand1(x,1,1);
%     
%     NeumIntegral1 = GaussQuad_Line3pt(x2,x1,integrand1x);
%    
%     GRHS(NIone(1)) = GRHS(NIone(1))-NeumIntegral1;
%     
%     %last element in top bndry elem 
%     NIend = TriMesh.topbndryelem(size(TriMesh.topbndryelem,1),:);
%     x1 = TriMesh.Xcoord(NIend(1));
%     y1 = TriMesh.Ycoord(NIend(1));
%     x2 = TriMesh.Xcoord(NIend(2));
%     y2 = TriMesh.Ycoord(NIend(2));
%     x3 = TriMesh.Xcoord(NIend(3));
%     y3 = TriMesh.Ycoord(NIend(3));
%     
%    localbasisfunct = localbasisfunct_CGP1(x1,y1,x2,y2,x3,y3);
% 
%     integrand2 = @(x,y,n)localbasisfunct.phi2(x,y).*Neumfunct(x,y,n);
%    
%     integrand2x = @(x)integrand2(x,1,1);
%     
%     NeumIntegral2 = GaussQuad_Line3pt(x2,x1,integrand2x);
%     
%     GRHS(NIend(2)) = GRHS(NIend(2))-NeumIntegral2;


%bottom boundary add Neumann BCs
for i =1:size(TriMesh.bottombndryelem,1)
    NI = TriMesh.bottombndryelem(i,:);
    x1 = TriMesh.Xcoord(NI(1));
    y1 = TriMesh.Ycoord(NI(1));
    x2 = TriMesh.Xcoord(NI(2));
    y2 = TriMesh.Ycoord(NI(2));
    x3 = TriMesh.Xcoord(NI(3));
    y3 = TriMesh.Ycoord(NI(3));
    
    areaT = 0.5*(x1*y1-x1*y3-y1*x2+y1*x3+x2*y3-y2*x3);
    Phi1 = @(x,y) (1/(2*areaT))*(y3*(x2-x)+y2*(x-x3)+y*(x3-x2));
    Phi2 = @(x,y) (1/(2*areaT))*(y3*(x-x1)+y1*(x3-x)+y*(x1-x3));
    
    %localbasisfunct_CGP1(x1,y1,x2,y2,x3,y3);
    
    integrand1 = @(x,y)Phi1(x,y).*Neumfunct_bottom(x,y);
    integrand2 = @(x,y)Phi2(x,y).*Neumfunct_bottom(x,y);
    
    %integrand1 = @(x,y)localbasisfunct.phi2(x,y).*Neumfunct_bottom(x,y);
    %integrand2 = @(x,y)localbasisfunct.phi1(x,y).*Neumfunct_bottom(x,y);
    
    integrand1x = @(x)integrand1(x,0);
    integrand2x = @(x)integrand2(x,0);
    
    NeumIntegral1 = GaussQuad_Line3pt(x1,x2,integrand1x);
    NeumIntegral2 = GaussQuad_Line3pt(x1,x2,integrand2x);

    GRHS(NI(1)) = GRHS(NI(1))-NeumIntegral1;
    GRHS(NI(2)) = GRHS(NI(2))-NeumIntegral2;
end

% %first elem in bottom boundary    AGAIN, THIS MAY NOT BE NECESSARY 
% NIone = TriMesh.bottombndryelem(1,:);
%     x1 = TriMesh.Xcoord(NIone(1));
%     y1 = TriMesh.Ycoord(NIone(1));
%     x2 = TriMesh.Xcoord(NIone(2));
%     y2 = TriMesh.Ycoord(NIone(2));
%     x3 = TriMesh.Xcoord(NIone(3));
%     y3 = TriMesh.Ycoord(NIone(3));
%     
%    localbasisfunct = localbasisfunct_CGP1(x1,y1,x2,y2,x3,y3);
%     
%     integrand2 = @(x,y,n)localbasisfunct.phi2(x,y).*Neumfunct(x,y,n);
%     
%     integrand2x = @(x)integrand2(x,0,-1);
%     
%     NeumIntegral2 = GaussQuad_Line3pt(x2,x1,integrand2x);
%    
%     GRHS(NIone(2)) = GRHS(NIone(2))-NeumIntegral2;
%     
%     %last element in top bndry elem 
%     NIend = TriMesh.topbndryelem(size(TriMesh.bottombndryelem,1),:);
%     x1 = TriMesh.Xcoord(NIend(1));
%     y1 = TriMesh.Ycoord(NIend(1));
%     x2 = TriMesh.Xcoord(NIend(2));
%     y2 = TriMesh.Ycoord(NIend(2));
%     x3 = TriMesh.Xcoord(NIend(3));
%     y3 = TriMesh.Ycoord(NIend(3));
%     
%    localbasisfunct = localbasisfunct_CGP1(x1,y1,x2,y2,x3,y3);
% 
%     integrand1 = @(x,y,n)localbasisfunct.phi1(x,y).*Neumfunct(x,y,n);
%    
%     integrand1x = @(x)integrand1(x,0,-1);
%     
%     NeumIntegral1 = GaussQuad_Line3pt(x2,x1,integrand1x);
%     
%     GRHS(NIend(1)) = GRHS(NIend(1))-NeumIntegral1;
% 
%     

end

