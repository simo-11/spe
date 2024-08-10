%% set parameters
rect_width=10;
rect_height=10;
rect_a=rect_width/2;
rect_b=rect_height/2;
rect_nMax=3;
rect_x0=rect_a;
rect_y0=rect_b;
%% create tps using increasing number of sites
figure(401)
for rand_count=3:30:500
    c_x = rect_width*[0 rand([1 rand_count]) 1];
    c_y = rect_height*[0 rand([1 rand_count]) 1];
    c_z = rect_psi(c_x,c_y,rect_x0,rect_y0,rect_nMax,rect_a,rect_b);
    tpaps_x=[c_x; c_y];
    tpaps_y=c_z;
    tps = tpaps(tpaps_x,tpaps_y,1);
    fnplt(tps);
    titletext=sprintf("TPS, site count=%d",rand_count+2);
    title(titletext);
    pause(0.2);
end
%% create fits using increasing number of points
models=["linearinterp" "cubicinterp" ...
    "biharmonicinterp" "thinplateinterp","lowess"];
for rand_count=3:30:500
    c_x = rect_width*[0 rand([1 rand_count]) 1];
    c_y = rect_height*[0 rand([1 rand_count]) 1];
    c_z = rect_psi(c_x,c_y,rect_x0,rect_y0,rect_nMax,rect_a,rect_b);
    for mi=1:models.size(2)
        figure(410+mi);
        [ff,gof,output] = fit([c_x',c_y'],c_z',models(mi));
        plot(ff);
        titletext=sprintf("%s, point count=%d",models(mi),rand_count+2);
        title(titletext);
        pause(0.1);
    end
end
%% plot warping function based on analytical solution
figure(403)
[X,Y]=meshgrid(0:rect_width/51:rect_width,...
    0:rect_height/51:rect_height);
Z=rect_psi(X,Y,rect_x0,rect_y0,rect_nMax,rect_a,rect_b);
surf(X,Y,Z);
title('Analytical solution')