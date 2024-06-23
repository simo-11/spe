function [z]=rect_psi(x1,y1,x0,y0,nMax,a,b)
% warping function for rectangle   
%{ 
en.wikiversity.org/wiki/Warping_functions
%}
if a<b
    error("a must be >= b")
end
x=x1-x0;
y=y1-y0;
z=x.*y;
sm=32*a^2/(pi^3);
for n=0:nMax
    kn=(2*n+1)*pi/(2*a);
    sv=(((-1)^n).*sin(kn.*x).*sinh(kn.*y))/...
        ((2*n+1)^3*cosh(kn*b));
    z=z-sm*sv;
end
end        
