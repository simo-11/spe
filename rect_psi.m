function [z]=rect_psi(x1,y1,x0,y0,nMax,a,b)
% warping function for rectangle
%{ 
en.wikiversity.org/wiki/Warping_functions
a=half width in x-direction
b=half height in y-direction
%}
x=x1-x0;
if x>a
    error("x (%f) cannot be larger than x1-x0 (%f-%f)",x,x1,x0);
end
y=y1-y0;
if y>b
    error("y (%f) cannot be larger than y1-y0 (%f-%f)",y,y1,y0);
end
%{
la=longer axis
sa=shorter axis
%}
if a>=b
    la=a;
    sa=b;
    sinf=x;
    sinhf=y;
else
    la=b;
    sa=a;
    sinf=y;
    sinhf=x;
end
z=x.*y;
sm=32*la^2/(pi^3);
for n=0:nMax
    kn=(2*n+1)*pi/(2*la);
    sv=(((-1)^n).*sin(kn.*sinf).*sinh(kn.*sinhf))/...
        ((2*n+1)^3*cosh(kn*sa));
    z=z-sm*sv;
end
end        
