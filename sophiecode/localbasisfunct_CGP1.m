function [ LocalBasisFuncts ] = localbasisfunct_CGP1( x1,y1,x2,y2,x3,y3 )
%a structure containing the local basis functions for one element in
%TriMesh

Tmtx = [1, x1, y1; 1, x2, y2; 1, x3, y3];
areaT = 0.5*(det(Tmtx));
c = 1/(2*areaT);

Phi1 = @(x,y) c*(y3*(x2-x)+y2*(x-x3)+y*(x3-x2));
Phi2 = @(x,y) c*(y3*(x-x1)+y1*(x3-x)+y*(x1-x3));
Phi3 = @(x,y) c*(y2*(x1-x)+y1*(x-x2)+y*(x2-x1));

LocalBasisFuncts = struct('phi1', Phi1, 'phi2', Phi2, 'phi3', Phi3);

end

