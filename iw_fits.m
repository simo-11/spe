function iw_fits(ao)
arguments
    ao.height=100
    ao.width=100
    ao.plot_fitted_surface=1
    ao.test_for_duplicates=1
    ao.save_xy_plot=1
    ao.n_max=4
end
rect_width=ao.width/1000;
rect_height=ao.height/1000;
test_for_duplicates=ao.test_for_duplicates;
plot_fitted_surfaces=ao.plot_fitted_surface;
save_xy_plot=ao.save_xy_plot;
rect_nMax=ao.n_max;
rect_a=rect_width/2;
rect_b=rect_height/2;
rect_x0=rect_a;
rect_y0=rect_b;
abs_tol=rect_width^3*rect_height^3/144;
point_set=haltonset(2,skip=30);
close('all')
figure(400)
[X,Y]=meshgrid(0:rect_width/51:rect_width,...
    0:rect_height/51:rect_height);
Z=rect_psi(X,Y,rect_x0,rect_y0,rect_nMax,rect_a,rect_b);
surf(X,Y,Z);
w=@(x,y)rect_psi(x,y,rect_x0,rect_y0,rect_nMax,rect_a,rect_b).^2;
Iw_a=integral2(w,0,rect_width,0,rect_height,...
      AbsTol=abs_tol,RelTol=0.0001);
titletext=sprintf("Analytical %Gx%G mm, n=%d, Iw=%.4G",...
    rect_width*1000,rect_height*1000,rect_nMax,Iw_a);
title(titletext);
da=daspect;
daspect([da(2) da(2) da(3)]);
models=["linearinterp" "cubicinterp" "poly44"...
     "thinplateinterp"];
% biharmonicinterp is similar as thinplateinterp
% naturalinterp is quite similar is linearinterp
line_specs=["--." "-o" "--x" "-^" "-v" ":o" ":x" "-."];
random_counts=[10:10:100 100:30:400];
edge_point_counts=5:2:9;% [12:3:30];
%x_values=4:1:17;% use to view shapes at low point counts
%x_values=4:17;
edge_size=size(edge_point_counts,2);
random_size=size(random_counts,2);
mp=get(0,'MonitorPositions');
fig_height=mp(4)/models.size(2)-20;
fig_width=mp(3)/5;
for ei=1:edge_size
    for mi=1:models.size(2)
        cmd=sprintf("%s_values=[];",models(mi));
        eval(cmd);
        cmd=sprintf("%s_x_values=[];",models(mi));
        eval(cmd);
    end
    edge_point_count=edge_point_counts(ei);
    d_e=1/(edge_point_count-1);
    e_p_0_1=linspace(0,1-d_e,edge_point_count-1);
    e_p_1_1=linspace(1,1,edge_point_count-1);
    e_p_1_0=linspace(1,0+d_e,edge_point_count-1);
    e_p_0_0=linspace(0,0,edge_point_count-1);
    for ri=1:random_size
        plot_line="plot(";
        rand_count=random_counts(ri);
        % create points on edges to help with some models
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
            ff = fit([c_x',c_y'],c_z',model);
            w=@(x,y)ff(x,y).^2;
            Iw=integral2(w,0,rect_width,0,rect_height,...
                AbsTol=abs_tol,RelTol=0.0001);
            ev=100*(Iw-Iw_a)/Iw_a; %#ok<NASGU> used using eval
            cmd=sprintf("%s_values(%d)=ev;",model,x_i);
            eval(cmd);
            if plot_fitted_surfaces
                fig=figure(410+mi);
                if ri==1
                    fig.Position=[0 (mi-1)*fig_height ...
                        fig_width fig_height];
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
        if ri==1
            fig.Position=[fig_width+edge_point_count*10 ...
                edge_point_count*50 ...
                fig_width fig_height];
        end
        plot_line=sprintf("%s)",plot_line);
        eval(plot_line);
        titletext=sprintf("I_w error with %d edge points",...
            edge_point_count);
        title(titletext);
        xlabel('total number of points');
        ylabel('Error %');
        ylim([-10 10]);
        %ylim("auto");
        yscale('linear');
        xscale('linear');
        legend(active_models(:));
        if ri==random_size && save_xy_plot
            filename=sprintf("Iw_error-%G-%G-%d",...
                rect_width*1000,rect_height*1000,edge_point_count);
            save_pdf_and_fig(fig,filename);
        end
        pause(0.1);
    end
end