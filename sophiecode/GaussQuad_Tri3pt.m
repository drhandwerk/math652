function [ GaussIntegral ] = GaussQuad_Tri3pt( x1,y1, x2,y2,x3,y3, integrandfunct )
%Gaussian quadrature on a triangle with three quadrature points
format long

weights = (1/3)*ones(3,1);
barycentric_coords = zeros(3,3);
barycentric_coords(1,1) = 2/3;
barycentric_coords(1,2) = 1/6;
barycentric_coords(1,3) = 1/6;

barycentric_coords(2,1) = 1/6;
barycentric_coords(2,2) = 2/3;
barycentric_coords(2,3) = 1/6;

barycentric_coords(3,1) = 1/6;
barycentric_coords(3,2) = 1/6;
barycentric_coords(3,3) = 2/3;

Tmtx = [1, x1, y1; 1, x2, y2; 1, x3, y3];
areaT = 0.5*(det(Tmtx));

quadpts = zeros(3,2);
for i = 1:3
    quadpts(i,1) = barycentric_coords(i,1)*x1+barycentric_coords(i,2)*x2+barycentric_coords(i,3)*x3;
    quadpts(i,2) = barycentric_coords(i,1)*y1+barycentric_coords(i,2)*y2+barycentric_coords(i,3)*y3;
end

GaussIntegral = 0;

for i = 1:3
    GaussIntegral = GaussIntegral + weights(i)*integrandfunct(quadpts(i,1),quadpts(i,2));
end

GaussIntegral = areaT*GaussIntegral;

end

