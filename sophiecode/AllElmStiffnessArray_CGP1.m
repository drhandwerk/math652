function [ AllESM, sumALLESM] = AllElmStiffnessArray_CGP1( TriMesh )
%input: a triangular mesh
%output: an array where each page corresponds to an element stiffness mtx

AllESM = zeros(3,3,TriMesh.numbelem);
sumALLESM = zeros(TriMesh.numbelem, 6);
for i = 1:TriMesh.numbelem
    %vector consisting of the nodes in element i
    NI = TriMesh.elem_node(i,:);  
    
    %Obtain the coordinates for the nodes on element i
    x1 = TriMesh.Xcoord(NI(1));
    y1 = TriMesh.Ycoord(NI(1));
    x2 = TriMesh.Xcoord(NI(2));
    y2 = TriMesh.Ycoord(NI(2));
    x3 = TriMesh.Xcoord(NI(3));
    y3 = TriMesh.Ycoord(NI(3));
 
    AllESM(:,:,i) = ElmStiffnessMtx_CGP1(x1,y1,x2,y2,x3,y3);
    sumc1 = sum(AllESM(:,1,i));
    sumc2 = sum(AllESM(:,2,i));
    sumc3 = sum(AllESM(:,3,i));
    sumr1 = sum(AllESM(1,:,i));
    sumr2 = sum(AllESM(2,:,i));
    sumr3 = sum(AllESM(3,:,i));
   
    sumALLESM(i,:)= [sumc1,sumc2,sumc3,sumr1, sumr2, sumr3];
end

    

end

