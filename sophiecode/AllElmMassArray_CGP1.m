function [ AllEMM ] = AllElmMassArray_CGP1( TriMesh )
%Creates an array of all of the element stiffness matricies for a
%trianglualr mesh. 
AllEMM = zeros(3,3,TriMesh.numbelem);

for i = 1:TriMesh.numbnodes
    NI = TriMesh.elem_node(i,:);
    x1 = TriMesh.Xcoord(NI(1));
    y1 = TriMesh.Ycoord(NI(1));
    x2 = TriMesh.Xcoord(NI(2));
    y2 = TriMesh.Ycoord(NI(2));
    x3 = TriMesh.Xcoord(NI(3));
    y3 = TriMesh.Ycoord(NI(3));
    
    AllEMM(:,:,i) = ElementMassMtx_CGP1(x1,y1,x2,y2,x3,y3);
end


end

