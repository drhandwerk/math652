function [ GSM ] = AssembleGSM_CGP1( TriMesh, AllElmStiffnessArray )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

GSM = zeros(TriMesh.numbnodes, TriMesh.numbnodes);

for i = 1:TriMesh.numbelem
    NI = TriMesh.elem_node(i,:);  %vector containing the counterclockwise ordering of node indicies for element i
    
    for j = 1:3
        for k = 1:3
            GSM(NI(j),NI(k)) = GSM(NI(j),NI(k)) + AllElmStiffnessArray(j,k,i);
        end
    end

end
%GSM
GSM = sparse(GSM);

end

