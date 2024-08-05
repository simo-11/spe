% verification of thin_plate_spline.py
%% create tps
rect_width=12;
rect_height=10;
rect_a=rect_width/2;
rect_b=rect_height/2;
rect_nMax=3;
rect_x0=rect_a;
rect_y0=rect_b;
c_x = [0, 1, 1, 2, 4];
c_y = [0, 1, 2, 2, 4];
c_z = rect_psi(c_x,c_y,rect_x0,rect_y0,rect_nMax,rect_a,rect_b);
tpaps_x=[c_x; c_y];
tpaps_y=c_z;
tps = tpaps(tpaps_x,tpaps_y,1);
%% Evaluate TPS at new points
x = linspace(0, 4, 5);
y = linspace(0, 4, 5);
z = fnval(tps,[x; y]);
%% plot warping function based on analytical solution
[X,Y]=meshgrid(0:rect_width/51:rect_width,...
    0:rect_height/51:rect_height);
Z=rect_psi(X,Y,rect_x0,rect_y0,rect_nMax,rect_a,rect_b);
surf(X,Y,Z);
