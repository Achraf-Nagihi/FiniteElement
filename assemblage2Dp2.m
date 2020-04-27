function [ A,F ] = assemblage2Dp2( alpha,beta,f,p,t,b )

P=p';
T=(t(1:6,:))';

E=size(T,1);%nombre d'éléments
N=size(P,1); %nombre de noeuds


%assemblage de la matrice
A=zeros(N,N);
F=zeros(N,1);
for i=1:E
    nodes=T(i,:);
    elemM=masse2dp2(i,p,t);
    elemK=raideur2dp2(i,p,t);
    elemA=alpha*elemM+beta*elemK;
    A(nodes,nodes)=A(nodes,nodes)+elemA;
    elemF=secondF2dp2(i,f,p,t);
    F(nodes)=F(nodes)+elemF ;
end

%les conditions aux limites
A(b,:)=0;
A(:,b)=0;
F(b)=0; 
A(b,b)=eye(length(b),length(b));

end

