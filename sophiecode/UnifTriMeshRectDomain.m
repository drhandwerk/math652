function [ TriMesh] = UnifTriMeshRectDomain( x1, x2, y1, y2, nx, ny )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

numbnodes = (nx+1)*(ny+1);
numbelem = 2*nx*ny;
numbedg = nx*(ny+1)+ny*(nx+1)+nx*ny;
numbvertedg   = ny*(nx+1);
numbhorizedg  = nx*(ny+1);
numbdiagedg = nx*ny;

hx = (x2-x1)/nx;
hy = (y2-y1)/ny;

x = linspace(x1, x2, nx+1);
y = linspace(y1, y2, ny+1);

%creates matricies with the bottom left node in location X(1,1) and Y(1,1).

[X,Y] = meshgrid(x,y);
meshgrid_mtx = X;
meshgrid_mtx(:,:,2) = Y;
%Node 1 in bottom left corner is X(1) and Y(1). Then node 2 is X(2),Y(2).
X = X(:);
Y = Y(:);
nodenumber = reshape((1:numbnodes),ny+1, nx+1);
nodes = [X,Y];

elem_node = zeros(numbelem, 3);

%saving elements vs nodes 
for i = 1: nx
    for j = 1:ny
    %odd numbered elements (left side of bisected rectangle) where
    %elem(ii,1) corresponds to the bottom left corner and proceed in
    %counterclockwise orientation
    ii = (i-1)*(2*ny)+2*j-1;
    elem_node(ii,1) = (i-1)*(ny+1)+j;
    elem_node(ii,2) = elem_node(ii,1)+(ny+1);
    elem_node(ii,3) = elem_node(ii,1)+1;
    
    %even numbered elements (right side of the bisected rectangle) where
    %elem(jj,1) corresponds to the top right corner and proceed in
    %counterclockwise orientation.
    jj = (i-1)*(2*ny)+2*j;
    elem_node(jj,1) = i*(ny+1)+j+1;
    elem_node(jj,2) = elem_node(jj,1)-(ny+1);
    elem_node(jj,3) = elem_node(jj,1)-1;
    
    end
end

edg_node = zeros(numbedg, 2);
%the first ny(nx+1) rows correspond to vertical edges, the middle
%nx(ny+1)rows correspond to the horizontal edges, and the last nx*ny rows
%correspond to the diagonal edges

%vertical edges
for i = 1: nx+1
    for j = 1: ny 
        k = (i-1)*ny + j;
        edg_node(k,1) = (i-1)*(ny+1)+j;
        edg_node(k,2) = edg_node(k,1)+1;
    end
end

%horizontal edges
for i = 1:nx 
    for j = 1: ny+1
        k = numbvertedg + (j-1)*nx+i;
        edg_node(k,1) = (i-1)*(ny+1)+j;
        edg_node(k,2) = edg_node(k,1) + (ny+1);
    end
end

%diag edges
for i = 1:nx
    for j =1:ny
        k = numbvertedg + numbhorizedg + (i-1)*ny+j;
        edg_node(k,1 ) = (i-1)*(ny+1)+j+1;
        edg_node(k,2 ) = edg_node(k,1) + ny;
    end
end      

%identify interior nodes and boundary nodes
intnodes = nodenumber(2:end-1, 2:end-1);
bndrynodes = setdiff(nodenumber(:),intnodes(:));

intnodes = intnodes(:);

%identify boundary elements
%left boundary elements
leftbndryelem = elem_node(1:2:2*ny,:);

rightbndryelem = elem_node(numbelem-2*ny+2:2:end, :);

topbndryelem = elem_node(2*ny:2*ny:end,:);

bottombndryelem = elem_node(1:2*ny:end,:);

TriMesh = struct('x1', x1, 'x2', x2, 'y1', y1, 'y2', y2, 'nx', nx, 'ny', ny, 'numbnodes', numbnodes, 'numbelem', numbelem, 'numbedg', numbedg,'numbvertedg', numbvertedg, 'numbhorizedg', numbhorizedg, 'numbdiagedg', numbdiagedg,'hx',hx, 'hy', hy, 'x', x, 'y', y, 'Xcoord', X, 'Ycoord',Y, 'nodes', nodes, 'nodeindex', nodenumber, 'elem_node', elem_node, 'edge_node', edg_node,'intnodes', intnodes,'bndrynodes', bndrynodes, 'leftbndryelem', leftbndryelem, 'rightbndryelem', rightbndryelem,'topbndryelem',topbndryelem, 'bottombndryelem',bottombndryelem, 'meshgridXY',meshgrid_mtx);

end

