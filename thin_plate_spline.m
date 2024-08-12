%% set parameters and plot warping function based on analytical solution
rect_width=0.1;
rect_height=0.1;
if rect_height > rect_width
    error("rect_width must be > rect_height")
end
rect_a=rect_width/2;
rect_b=rect_height/2;
rect_nMax=3;
rect_x0=rect_a;
rect_y0=rect_b;
abs_tol=rect_width^3*rect_height^3/144;
include_lowess=0; % it is not suitable for this case 
point_set=haltonset(2,skip=1);
figure(400)
[X,Y]=meshgrid(0:rect_width/51:rect_width,...
    0:rect_height/51:rect_height);
Z=rect_psi(X,Y,rect_x0,rect_y0,rect_nMax,rect_a,rect_b);
surf(X,Y,Z);
w=@(x,y)rect_psi(x,y,rect_x0,rect_y0,rect_nMax,rect_a,rect_b).^2;
Iw_a=integral2(w,0,rect_width,0,rect_height,...
      AbsTol=abs_tol,RelTol=0.0001);
titletext=sprintf("Analytical Iw=%.4G",Iw_a);
title(titletext);
%% create tps using increasing number of sites
figure(401)
for rand_count=3:30:500
    h_set=net(point_set,rand_count);
    c_x = rect_width* [0 1 1 0 h_set(:,1)'];
    c_y = rect_height*[0 0 0 1 h_set(:,2)'];
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
    "biharmonicinterp" "thinplateinterp"];
if include_lowess
    models=[models "lowess"]; %#ok<UNRCH>
end
x_values=[20:5:100 100:30:501];
x_len=size(x_values,2);
for mi=1:models.size(2)
    cmd=sprintf("%s_values=zeros(1,%d);",models(mi),x_len);
    eval(cmd);
end
plot_line="plot(";
for ri=1:x_len
    rand_count=x_values(ri);
    % create points on edges to help with some models
    edge_point_count=floor(sqrt(rand_count));
    d_e=1/(edge_point_count-1);
    e_p_0_1=linspace(0,1-d_e,edge_point_count-1);
    e_p_1_1=linspace(1,1,edge_point_count);
    e_p_1_0=linspace(1,0+d_e,edge_point_count-1);
    e_p_0_0=linspace(0,0,edge_point_count);
    h_set=net(point_set,rand_count);
    c_x = rect_width* [e_p_0_1 e_p_1_1 e_p_1_0 e_p_0_0 h_set(:,1)'];
    c_y = rect_height*[e_p_0_0 e_p_0_1 e_p_1_1 e_p_1_0 h_set(:,2)'];
    c_z = rect_psi(c_x,c_y,rect_x0,rect_y0,rect_nMax,rect_a,rect_b);
    for mi=1:models.size(2)
        figure(410+mi);
        [ff,gof,output] = fit([c_x',c_y'],c_z',models(mi));
        plot(ff);
        da=daspect;
        daspect([da(2) da(2) da(3)]);
        w=@(x,y)ff(x,y).^2;
        Iw=integral2(w,0,rect_width,0,rect_height,...
            AbsTol=abs_tol,RelTol=0.0001);
        ev=100*abs((Iw-Iw_a)/Iw_a);
        cmd=sprintf("%s_values(%d)=ev;",models(mi),ri);
        eval(cmd);
        titletext=sprintf("%s, point count=%d, Iw=%.4G",...
            models(mi),size(c_x,2),Iw);
        title(titletext);
        pause(0.1);
        if rand_count==x_values(1)
            if mi==1
                comma="";
            else
                comma=",";
            end
            plot_line=sprintf("%s%sx_values,%s_values",...
                plot_line,comma,models(mi));
        end
    end
    figure(420);
    if rand_count==x_values(1)
        plot_line=sprintf("%s)",plot_line);
    end
    eval(plot_line);
    title('I_w error');
    xlabel('number of points');
    ylabel('Error %');
    yscale('log');
    xscale('log');
    legend(models(:));
    pause(0.1);
end
