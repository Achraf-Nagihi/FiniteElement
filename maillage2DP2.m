function [ p2,t2,b2 ] = maillage2DP2( p,t,e )

%génération du maillage des milieux
t=(t(1:3,:))';
T=size(t,1);
t0=zeros(T,3);
t2=[t,t0];
nmbp=length(p(1,:));
p0=zeros(2,1);
p2=[p,p0];
k=nmbp+1;
for i=1:T
    nod=t2(i,1:3)'; %les noeuds associé à l'élement
    %définition des coordonées des noeuds "nod"
    x1=p(1,nod(1));   y1=p(2,nod(1));
    x2=p(1,nod(2));   y2=p(2,nod(2));
    x3=p(1,nod(3));   y3=p(2,nod(3));
    %calcul des milieux
    mx1=(x1+x2)/2;    my1=(y1+y2)/2;
    mx2=(x2+x3)/2;    my2=(y2+y3)/2;
    mx3=(x3+x1)/2;    my3=(y3+y1)/2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dans cette partie on va remplir le tableau des noeuds au milieu tout%
    % assurant qu'il n'esxiste pas de redondance                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    test1=0;
    n=length(p2(1,:));
    for j=(nmbp+1):n
        if p2(1,j)==mx1 && p2(2,j)==my1
            t2(i,4)=j;
            k=k-1;
            test1=1;
        end
    end
        if test1==0;
            p2(1,k)=mx1;
            p2(2,k)=my1;
            t2(i,4)=k;
        end
    %%%%%%%%%%%%%%%%%%%
    test2=0;
    n=length(p2(1,:));
    for j=(nmbp+1):n
        if p2(1,j)==mx2 && p2(2,j)==my2
            t2(i,5)=j;
            k=k-1;
            test2=1;
        end
    end
        if test2==0
            p2(1,k+1)=mx2;
            p2(2,k+1)=my2;
            t2(i,5)=k+1;
        end
    %%%%%%%%%%%%%%%%
    test3=0;
    n=length(p2(1,:));
    for j=(nmbp+1):n
        if p2(1,j)==mx3 && p2(2,j)==my3
            t2(i,6)=j;
            k=k-1;
            test3=1;
        end 
    end
        if test3==0
            p2(1,k+2)=mx3;
            p2(2,k+2)=my3;
            t2(i,6)=k+2;
        end
        k=k+3;
end

%les éléments du bord
k=1;
bm(1,1)=0;
a=(p(1,e(1,:))+p(1,e(2,:)))/2;
b=(p(2,e(1,:))+p(2,e(2,:)))/2;
for j=1:length(a)
    for i=1:length(p2(1,:))
        if a(j)==p2(1,i) && b(j)==p2(2,i)
            bm(1,k)=i;
            k=k+1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%résultats demandés
t2=t2';
b2=[e(1,:) bm];
%p2 est déja prêt
end


