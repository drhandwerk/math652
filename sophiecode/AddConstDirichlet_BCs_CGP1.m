function [ GSM_DBC, GRHS_DBC ] = AddConstDirichlet_BCs_CGP1( GRHS, GSM, TriMesh, DirichletConst )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[ leftandright, top , bottom] = Dirichlet_Neumann_BCs_rectdomain_CGP1_hw2( TriMesh );
%leftandright corresponds to the dirichlet bndry nodes for m652 hw2
topbottom = vertcat(top',bottom');

for i = 1: size(leftandright,1)
    diag = GSM(leftandright(i),leftandright(i));
    GSM(leftandright(i),:)=0;
    GSM(leftandright(i),leftandright(i))=diag;
    GRHS(leftandright(i))=DirichletConst;
end

%this part is to include dirichlet BCs on whole domain
for i = 1: size(topbottom,1)
    diag = GSM(topbottom(i),topbottom(i));
    GSM(topbottom(i),:)=0;
    GSM(topbottom(i),topbottom(i))=diag;
    GRHS(topbottom(i))=0;
end

GSM_DBC = GSM;

GRHS_DBC = GRHS;

end

