%% Make a mesh over L-shaped domain with initial edges size h0. 
% Save the nodes, triangles, and boundary nodes. 
function [p,t,b] = genLmesh(h0)
addpath('../../distmesh')
fd = @(p) ddiff(drectangle(p,-1,1,-1,1),drectangle(p,0,1,-1,0));
[p,t] = distmesh2d(fd,@huniform,h0,[-1,-1;1,1],[0,0;1,0;1,1;-1,1;-1,-1;0,-1]);
e = boundedges(p,t);
b = union(e(:,1),e(:,2));

save('distmeshdata.mat','p','t','b')
end

