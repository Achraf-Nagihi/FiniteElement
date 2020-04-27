function [ p,t,e ] = maillage2DP1( x1,y1,x2,y2,x3,y3,x4,y4,h )
rect = [3 ;4 ;x1; x2; x3; x4; y1; y2; y3; y4];
dl = decsg(rect);
[p,e,t] = initmesh(dl,'hmax',h);
end

