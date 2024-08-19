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
test_for_duplicates=1;
plot_fitted_surfaces=0;
include_lowess=0; % lowess is not suitable for this case 
point_set=haltonset(2,skip=30);
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
da=daspect;
daspect([da(2) da(2) da(3)]);
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
models=["linearinterp" "cubicinterp" "poly44"...
    "biharmonicinterp" "thinplateinterp"];
line_specs=["--." "-o" "--x" "-^" "-v" ":o" ":x" "-."];
if include_lowess
    models=[models "lowess"]; %#ok<UNRCH>
end
x_values=[20:5:100 100:30:501]; %#ok<NASGU>
%x_values=4:1:17;% use to view shapes at low point counts
%x_values=4:17;
x_len=size(x_values,2);
for mi=1:models.size(2)
    cmd=sprintf("%s_values=[];",models(mi));
    eval(cmd);
    cmd=sprintf("%s_x_values=[];",models(mi));
    eval(cmd);
end
mp=get(0,'MonitorPositions');
fig_height=mp(4)/models.size(2)-20;
fig_width=mp(3)/5;
px_values=[];
for ri=1:x_len
    plot_line="plot(";
    rand_count=x_values(ri);
    % create points on edges to help with some models
    edge_point_count=floor(sqrt(rand_count));
    %edge_point_count=5;
    d_e=1/(edge_point_count-1);
    e_p_0_1=linspace(0,1-d_e,edge_point_count-1);
    e_p_1_1=linspace(1,1,edge_point_count-1);
    e_p_1_0=linspace(1,0+d_e,edge_point_count-1);
    e_p_0_0=linspace(0,0,edge_point_count-1);
    h_set=net(point_set,rand_count);
    c_x = rect_width* [e_p_0_1 e_p_1_1 e_p_1_0 e_p_0_0 h_set(:,1)'];
    c_y = rect_height*[e_p_0_0 e_p_0_1 e_p_1_1 e_p_1_0 h_set(:,2)'];
    if test_for_duplicates
        A=[c_x;c_y]';
        [B,BG]=groupcounts(A);
        if any(B>1)
            fprintf("Duplicates found\n");
            BM=cell2mat(BG);
            duplicate_rows=find(B>1);
            for di=duplicate_rows'
                fprintf("count=%d x=%g, y=%g\n",...
                    B(di),BM(di,1),BM(di,2));
            end    
            error("x,y duplicates found, see details above")
        end
    end
    c_z = rect_psi(c_x,c_y,rect_x0,rect_y0,rect_nMax,rect_a,rect_b);
    active_models=strings(0);
    for mi=1:models.size(2)
        model=models(mi);
        aFittype=fittype(model);
        if numcoeffs(aFittype)>size(c_x,2)
            continue;
        end
        active_models(size(active_models,2)+1)=model;
        cmd=sprintf("x_i=size(%s_x_values,2)+1;",model);
        eval(cmd);
        cmd=sprintf("%s_x_values(%d)=size(c_x,2);",model,x_i);
        eval(cmd);
        [ff,gof,output] = fit([c_x',c_y'],c_z',model);
        w=@(x,y)ff(x,y).^2;
        Iw=integral2(w,0,rect_width,0,rect_height,...
            AbsTol=abs_tol,RelTol=0.0001);
        ev=100*(Iw-Iw_a)/Iw_a;
        cmd=sprintf("%s_values(%d)=ev;",model,x_i);
        eval(cmd);
        if plot_fitted_surfaces
            fig=figure(410+mi); %#ok<UNRCH>% if not plotting
            if ri==1
                fig.Position=[0 (mi-1)*fig_height fig_width fig_height];
                fig.MenuBar='figure';
                fig.DockControls="on";
            end
            plot(ff,[c_x',c_y'],c_z');
            da=daspect;
            daspect([da(2) da(2) da(3)]);
            titletext=sprintf("%s, point count=%d, Iw=%.4G",...
                model,size(c_x,2),Iw);
            title(titletext);
            pause(0.1);
        end
        if mi==1
            comma="";
        else
            comma=",";
        end
        plot_line=sprintf('%s%s%s_x_values,%s_values,"%s"',...
            plot_line,comma,model,model,line_specs(mi));
    end
    fig=figure(420+edge_point_count);
    fig.Position=[fig_width edge_point_count*20 ...
        fig_width fig_height];
    plot_line=sprintf("%s)",plot_line);
    eval(plot_line);
    titletext=sprintf("I_w error with %d edge points",...
        edge_point_count);
    title(titletext);
    xlabel('number of points');
    ylabel('Error %');
    %ylim([-25 10]);
    ylim("auto");
    yscale('linear');
    xscale('linear');
    legend(active_models(:));
    pause(0.1);
end
