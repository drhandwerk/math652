function [ GaussLineInt ] = GaussQuad_Line3pt( x1, x2, integrandfunct )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
format long
weights = [5/18, 8/18, 5/18];

barycentric_coords = zeros(3,2);
barycentric_coords(1,1) = (1/2)*(1-sqrt(3/5));
barycentric_coords(1,2) = (1/2)*(1+sqrt(3/5));
barycentric_coords(2,1) = 1/2;
barycentric_coords(2,2) = 1/2;
barycentric_coords(3,1) = (1/2)*(1+sqrt(3/5));
barycentric_coords(3,2) = (1/2)*(1-sqrt(3/5));

linelength = x2-x1;

GaussLineInt = 0;
for i = 1:3
    GaussLineInt = GaussLineInt + weights(i)*integrandfunct(barycentric_coords(i,1)*x1+barycentric_coords(i,2)*x2);
end
GaussLineInt = linelength*GaussLineInt;



end

