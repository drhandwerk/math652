function [ localbasisgrad ] = localbasisgradfunct_CGP1( x1,y1,x2,y2,x3,y3 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Tmtx = [1, x1, y1; 1, x2, y2; 1, x3, y3];
areaT = 0.5*(det(Tmtx));

c = (1/(2*areaT));

gradphi1 = c*[y2-y3, x3-x2]';
gradphi2 = c*[y3-y1, x1-x3]';
gradphi3 = c*[y1-y2, x2-x1]';

localbasisgrad = struct('gradphi1', gradphi1,'gradphi2', gradphi2, 'gradphi3', gradphi3);

end

