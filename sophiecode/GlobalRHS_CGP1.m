function [ GRHS ] = GlobalRHS_CGP1( localRHS_CGP1, TriMesh )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


GRHS = zeros(TriMesh.numbnodes,1);

for i = 1: TriMesh.numbelem
    NI = TriMesh.elem_node(i,:);
    for j = 1:3
        GRHS(NI(j))=GRHS(NI(j)) + localRHS_CGP1(1,j,i);
    end
end

end

