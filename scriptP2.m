  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  Programme de résolution d'equation différentielle aux dérivées  %
  %  partielles    -alpha*u"+beta*u=f       sur le domaine           %
  %                       u=0               condition au bord        %
  %                                                                  %
  %      Réalisé par : MOUTATAHIR younes et TRADY jamal              %
  %                                                                  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



g=@(x,y) (x+1.5)*x*y*(y-1);  %fonction de test exacte
f=@(x,y) (x^2)*(y^2)-(x^2)*y+(1.5*x*(y^2))-(1.5*x*y)-(2*x^2)-(2*y^2)+(2*y)-(3*x);   %second memebre pour alpha=1 et beta=1
%f=@(x,y) -(2*x^2)-(2*y^2)+(2*y)-(3*x);      % second membre pour alpha=1 et beta=0
alpha=1;
beta=1;

x1=-1.5;y1=0;x2=0;y2=0;x3=0;y3=1;x4=-1.5;y4=1;h=0.06;
[ p,t,e ] = maillage2DP1( x1,y1,x2,y2,x3,y3,x4,y4,h );

[ p2,t2,b2 ] = maillage2DP2( p,t,e );
[ A,F ] = assemblage2Dp2( alpha,beta,f,p2,t2,b2 );

U=A\F    %le vecteur de la solution approchée
x=p2(1,:);  % les coodonées x
y=p2(2,:);% les coordonées y
UE=arrayfun(g,x,y);   %le vecteur solution exacte
ER=U-UE';   % le vecteur erreur pour verifier l'exactitude du programme

% le dessin de la solution
figure(1);
tri = delaunay(x,y);
trisurf(tri,x,y,U);
title('Representation de la solution approchee en P2')

%on peut aussi dessiner la fonction exacte pour la comparer avec l'approchée

figure(2);
trisurf(tri,x,y,UE);
title('Representation de la solution exacte en P2')


%dessin de l'erreur pour la comparer avec le plan 0
figure(3);
trisurf(tri,x,y,ER);
title('Representation de l erreur en P2')
colorbar;
view(2);