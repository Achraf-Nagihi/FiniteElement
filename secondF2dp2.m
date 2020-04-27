function [ elemF ] = secondF2dp2( i,f,p,t )


%définition des noeuds associé à l'élément :

x1=p(1,t(1,i));
x2=p(1,t(2,i));
x3=p(1,t(3,i));
y1=p(2,t(1,i));
y2=p(2,t(2,i));
y3=p(2,t(3,i));

%le jacobien :

J=[x2-x1,x3-x1;y2-y1,y3-y1];


%les paramètres de la quadrature de gauss 2d :

a1=1/3;
b1=1/3;
a2=0.2;
b2=0.6;
a3=0.2;
b3=0.2;
a4=0.6;
b4=0.2;
w(1)=-27/96;
w(2)=25/96;
w(3)=25/96;
w(4)=25/96;

%les fonctions de formes 

ph1=@(s,t) 2*s^2+2*t^2+4*s*t-3*s-3*t+1;
ph2=@(s,t) 2*s^2-s;
ph3=@(s,t) 2*t^2-t;
ph4=@(s,t) -4*s^2-4*s*t+4*s;
ph5=@(s,t) 4*s*t;
ph6=@(s,t) -4*t^2-4*s*t+4*t;

R1=[x1,y1]+(J*[a1,b1]')';
R2=[x1,y1]+(J*[a2,b2]')';
R3=[x1,y1]+(J*[a3,b3]')';
R4=[x1,y1]+(J*[a4,b4]')';

% calcule de l'intégrale par la quadrature de gauss

elemF(1,1)=(w(1)*f(R1(1,1),R1(1,2))*ph1(a1,b1))+(w(2)*f(R2(1,1),R2(1,2))*ph1(a2,b2))+(w(3)*f(R3(1,1),R3(1,2))*ph1(a3,b3))+(w(4)*f(R4(1,1),R4(1,2))*ph1(a4,b4));
elemF(2,1)=(w(1)*f(R1(1,1),R1(1,2))*ph2(a1,b1))+(w(2)*f(R2(1,1),R2(1,2))*ph2(a2,b2))+(w(3)*f(R3(1,1),R3(1,2))*ph2(a3,b3))+(w(4)*f(R4(1,1),R4(1,2))*ph2(a4,b4));
elemF(3,1)=(w(1)*f(R1(1,1),R1(1,2))*ph3(a1,b1))+(w(2)*f(R2(1,1),R2(1,2))*ph3(a2,b2))+(w(3)*f(R3(1,1),R3(1,2))*ph3(a3,b3))+(w(4)*f(R4(1,1),R4(1,2))*ph3(a4,b4));
elemF(4,1)=(w(1)*f(R1(1,1),R1(1,2))*ph4(a1,b1))+(w(2)*f(R2(1,1),R2(1,2))*ph4(a2,b2))+(w(3)*f(R3(1,1),R3(1,2))*ph4(a3,b3))+(w(4)*f(R4(1,1),R4(1,2))*ph4(a4,b4));
elemF(5,1)=(w(1)*f(R1(1,1),R1(1,2))*ph5(a1,b1))+(w(2)*f(R2(1,1),R2(1,2))*ph5(a2,b2))+(w(3)*f(R3(1,1),R3(1,2))*ph5(a3,b3))+(w(4)*f(R4(1,1),R4(1,2))*ph5(a4,b4));
elemF(6,1)=(w(1)*f(R1(1,1),R1(1,2))*ph6(a1,b1))+(w(2)*f(R2(1,1),R2(1,2))*ph6(a2,b2))+(w(3)*f(R3(1,1),R3(1,2))*ph6(a3,b3))+(w(4)*f(R4(1,1),R4(1,2))*ph6(a4,b4));

elemF=det(J)*elemF;

end