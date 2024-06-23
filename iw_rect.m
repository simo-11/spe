function v=iw_rect(ao)
arguments
    ao.height=0.1
    ao.width=0.1
    ao.debugLevel=0
    ao.n=3
end
if ao.width<ao.height
    error("width must be >= height")
end
a=ao.width/2;
b=ao.height/2;
x0=0;
y0=0;
w=@(x,y)rect_psi(x,y,x0,y0,ao.n,a,b).^2;
v=integral2(w,-a,a,-b,b);
end
