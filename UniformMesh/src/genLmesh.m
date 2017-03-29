function [p,t] = genLmesh()
addpath('../../distmesh')
fd = @(p) ddiff(drectangle(p,-1,1,-1,1),drectangle(p,0,1,-1,0));
[p,t]=distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[0,0;1,0;1,1;-1,1;-1,-1;0,-1]);
end

