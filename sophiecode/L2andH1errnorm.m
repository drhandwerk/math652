function [ L2norm, H1norm , H1seminorm] = L2andH1errnorm( TriMesh, sln, exactsln_funct, ExactSln_grad)
%inputs::
%TriMesh: triangular mesh on rectangular domain
%sln: vector of u-values on the nodes of TriMesh
%exactsln_funct and exactsln_grad: functions that take in x and y values
%output: the L2 norm of the error exactsln - FE sln  and the H1 norm of the
%error 

L2norm=0;
for i = 1: TriMesh.numbelem
    
    NI = TriMesh.elem_node(i,:);
    x1 = TriMesh.Xcoord(NI(1));
    y1 = TriMesh.Ycoord(NI(1));
    x2 = TriMesh.Xcoord(NI(2));
    y2 = TriMesh.Ycoord(NI(2));
    x3 = TriMesh.Xcoord(NI(3));
    y3 = TriMesh.Ycoord(NI(3));
    
    localbasisfunct = localbasisfunct_CGP1(x1,y1,x2,y2,x3,y3);
    
    u_h = @(x,y) (1)*(sln(NI(1))*localbasisfunct.phi1(x,y) + sln(NI(2))*localbasisfunct.phi2(x,y) + sln(NI(3))*localbasisfunct.phi3(x,y));
    
    integrand = @(x,y) (abs(exactsln_funct(x,y) -u_h(x,y))).^2;
    
    err_tri = GaussQuad_Tri3pt(x1,y1,x2,y2,x3,y3, integrand);
    
    L2norm = L2norm + err_tri;
end

 L2norm = sqrt(L2norm);

 H1semi=0;
 
for j = 1 : TriMesh.numbelem
    NI = TriMesh.elem_node(j,:);
    x1 = TriMesh.Xcoord(NI(1));
    y1 = TriMesh.Ycoord(NI(1));
    x2 = TriMesh.Xcoord(NI(2));
    y2 = TriMesh.Ycoord(NI(2));
    x3 = TriMesh.Xcoord(NI(3));
    y3 = TriMesh.Ycoord(NI(3));
    
    localbasisgrad = localbasisgradfunct_CGP1(x1,y1,x2,y2,x3,y3);
    
    grad_u_h = sln(NI(1)).*localbasisgrad.gradphi1 + sln(NI(2)).*localbasisgrad.gradphi2 + sln(NI(3)).*localbasisgrad.gradphi3;
    
    integrand_grad = @(x,y) (sqrt((ExactSln_grad.partialx(x,y) - grad_u_h(1))^2 + (ExactSln_grad.partialy(x,y)-grad_u_h(2))^2))^2;
    
    H1semi_tri = GaussQuad_Tri3pt(x1,y1,x2,y2,x3,y3, integrand_grad);
    
    H1semi = H1semi +H1semi_tri;
    
    H1seminorm = sqrt(H1semi);
    H1norm = H1semi + (L2norm)^2;
end

   H1norm = sqrt(H1norm);


end

