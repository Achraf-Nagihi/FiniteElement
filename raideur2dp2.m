function [ elemK ] = raideur2dp2(i,p,t)


%définition des noeuds associé à l'élément:

x1=p(1,t(1,i));
x2=p(1,t(2,i));
x3=p(1,t(3,i));
y1=p(2,t(1,i));
y2=p(2,t(2,i));
y3=p(2,t(3,i));

%le jacobien :

J=[x2-x1,x3-x1;y2-y1,y3-y1];

%les fonctions de formes :

ph1=@(s,t) 2*s^2+2*t^2+4*s*t-3*s-3*t+1;
ph2=@(s,t) 2*s^2-s;
ph3=@(s,t) 2*t^2-t;
ph4=@(s,t) -4*s^2-4*s*t+4*s;
ph5=@(s,t) 4*s*t;
ph6=@(s,t) -4*t^2-4*s*t+4*t;


%Calcul de l'integrale par la quadrature de Gauss 

elemK=det(J)*[gauss2dp2(ph1,ph1),gauss2dp2(ph1,ph2),gauss2dp2(ph1,ph3),gauss2dp2(ph1,ph4),gauss2dp2(ph1,ph5),gauss2dp2(ph1,ph6);
              gauss2dp2(ph2,ph1),gauss2dp2(ph2,ph2),gauss2dp2(ph2,ph3),gauss2dp2(ph2,ph4),gauss2dp2(ph2,ph5),gauss2dp2(ph2,ph6);
              gauss2dp2(ph3,ph1),gauss2dp2(ph3,ph2),gauss2dp2(ph3,ph3),gauss2dp2(ph3,ph4),gauss2dp2(ph3,ph5),gauss2dp2(ph3,ph6);
              gauss2dp2(ph4,ph1),gauss2dp2(ph4,ph2),gauss2dp2(ph4,ph3),gauss2dp2(ph4,ph4),gauss2dp2(ph4,ph5),gauss2dp2(ph4,ph6);
              gauss2dp2(ph5,ph1),gauss2dp2(ph5,ph2),gauss2dp2(ph5,ph3),gauss2dp2(ph5,ph4),gauss2dp2(ph5,ph5),gauss2dp2(ph5,ph6);
              gauss2dp2(ph6,ph1),gauss2dp2(ph6,ph2),gauss2dp2(ph6,ph3),gauss2dp2(ph6,ph4),gauss2dp2(ph6,ph5),gauss2dp2(ph6,ph6)];
          

end
