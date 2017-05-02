function [ ESM ] = ElmStiffnessMtx_CGP1( x1, y1, x2, y2, x3, y3 )
%inputs: corners of a triangle (x1,y1), (x2, y2), (x3, y3) (listed in
%counterclockwise ordering)

%compute area of the entire triagular element
Tmtx = [1, x1, y1; 1, x2, y2; 1, x3, y3];
areaT = 0.5*(det(Tmtx));
%areaT = 0.5*(x1*y1-x1*y3-y1*x2+y1*x3+x2*y3-y2*x3);
c = (1/(4*areaT));

ESM = zeros(3,3);

ESM(1,1) = c*((y2-y3)^2+(x3-x2)^2);
ESM(2,2) = c*((y3-y1)^2+(x1-x3)^2);
ESM(3,3) = c*((y1-y2)^2+(x2-x1)^2);

ESM(1,2) = c*((y2-y3)*(y3-y1)+(x3-x2)*(x1-x3));
ESM(2,1) = ESM(1,2);

ESM(1,3) = c*((y1-y2)*(y2-y3)+(x3-x2)*(x2-x1));
ESM(3,1) = ESM(1,3);

ESM(2,3) = c*((y3-y1)*(y1-y2)+(x1-x3)*(x2-x1));
ESM(3,2) = ESM(2,3);

end

