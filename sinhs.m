function z=sinhs(x0,y0,c1,c3,c5,xe,ye,H)
x=xe-x0;
y=ye-y0;
t0=x.*y;
t1=c1*sin(pi*y/H).*sinh(pi*x/H);
t3=c3*sin(3*pi*y/H).*sinh(3*pi*x/H);
t5=c5*sin(5*pi*y/H).*sinh(5*pi*x/H);
z=t0+t1+t3+t5;
end