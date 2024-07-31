function [c]=testRHS(ao)
arguments
    ao.height=150
    ao.width=150
    ao.t=8
    ao.r=16 % outer radius
    ao.n_r=8
    ao.models=["cubicinterp","tps"]
    ao.cubs=["rbfcub"]
    ao.scat_type='halton'
    ao.cards=[10,20,25,30]    
    ao.debugLevel=0
    ao.max_epc=100
    ao.plot=0
    % bitwise and
    % 1:domain 
    % 2:integration points 
    % 4:fitted surface
    ao.rsquareMin=0.2
    ao.check_area=1
    ao.latex=1
    ao.max_area_error_percent=0.1
end
%{
Test warping function fit using csv files in gen directory
%}
fp=sprintf("gen/warping-rhs-%g-%g-%g-%g-%g-*.csv",...
    ao.height,ao.width,ao.t,ao.r,ao.n_r);
list=dir(fp);
n=size(list,1);
c=cell(n,1);
ms=size(ao.models,2);
cs=size(ao.cubs,2);
H=ao.height/1000;
W=ao.width/1000;
T=ao.t/1000;
% define one quarter
if ao.r==0
    XV=[0 W/2 W/2 T T   0  ];
    YV=[0 0   T   T H/2 H/2];
    ps1=polyshape(XV,YV);
else
    [XV,YV]=getVertices(ao.r/1000-T,ao.n_r,H/2,W/2,T);
    ps1=polyshape(XV,YV);
end
ps2=rotate(ps1,90);
ps3=translate(ps2,[W 0]);
ps4=union(ps1,ps3);
ps5=rotate(ps4,180);
ps6=translate(ps5,[W H]);
ps=union(ps4,ps6);
domain.vertices=ps.Vertices;
domain.polyshape=ps;
domain.domain='rectangle';
[xlimit,ylimit]=boundingbox(domain.polyshape);
domain.dbox=[xlimit; ylimit];
cao.debugLevel=ao.debugLevel;
cao.plot=bitand(ao.plot,2);
if bitand(ao.plot,1)
    fig=gcf;
    if isempty(fig.Name)
        fig.Name=fp;
    else
        pos=fig.Position;
        fig=figure;
        fig.Name=fp;
        fig.Position(1)=pos(1)+pos(3);
    end
    plot(domain.polyshape);
    axis equal;
end
for i=1:n
    fn=list(i).name;
    fprintf("file=%s\n",fn);
    file=sprintf("%s/%s",list(i).folder,fn);
    t=readtable(file);
    rfn=replace(fn,"warping","results");
    rfn=replace(rfn,"csv","json");
    fn=sprintf("%s/%s",list(i).folder,rfn);
    cao.spr=jsondecode(fileread(fn));
    if ao.debugLevel>1
        fprintf("x values: %.3g - %.3g\n",min(t.x),max(t.x));
        fprintf("y values: %.3g - %.3g\n",min(t.y),max(t.y));
        fprintf("w values: %.3g - %.3g\n",min(t.w),max(t.w));
    end
    maxIw=(max(t.x)-min(t.x))*(max(t.y)-min(t.y))*(max(t.w)-min(t.w))^2;
    o.t=t;
    o.file=file;
    for mi=1:ms
        model=ao.models(mi);
        fitMethod='fit';
        if startsWith(model,"tps")
            fitMethod='tpaps';
        else  
            ft=model;
        end
        tic;
        cpu_start=cputime; %#ok<NASGU>
        switch fitMethod
            case 'fit'
                [f,gof,output,warnstr,errstr,convmsg]=...
                    fit([t.x t.y],t.w,ft); %#ok<ASGLU>
                if (gof.rsquare<ao.rsquareMin)
                    fprintf(['Model %s for file %s rejected,'...
                    ' rsquare=%.3g, min=%.3g\n'],...
                    model,fn,...
                    gof.rsquare,ao.rsquareMin);
                    continue;
                end                
                w=@(x,y)f(x,y);
            case 'tpaps'
                if strlength(model)>3
                    pin=str2double(extractAfter(model,3));
                else
                    pin=1;
                end
                f=tpaps([t.x t.y]',t.w',pin);
                w=@(x,y)reshape(fnval(f,[x(:)';y(:)']),...
                    size(x,1),[]);
        end
        es=sprintf("o.%s_fit_walltime=toc;",model);
        eval(es);
        es=sprintf("o.%s_fit_cputime=cputime-cpu_start;",model);
        eval(es);
        if ao.debugLevel>1
            disp(f);
            switch fitMethod
                case 'fit'
                disp(gof);
            end
        end
        for ci=1:cs
            cub=ao.cubs(ci);
            for card=ao.cards
                tic
                cao.card=card;
                [cao.centers,cao.dbox,cao.area_domain]=...
                   define_scattered_pointset(card,domain,ao.scat_type);            
                w_at_centers=w(cao.centers(:,1),cao.centers(:,2));
                weights=get_weights(domain,cub,cao);
                A=cao.spr.area;
                if ao.check_area
                    % Check domain definition
                    Ac=sum(weights);
                    error_percent=100*abs((Ac-A)/A);
                    if error_percent>ao.max_area_error_percent
                        fprintf(['Domain definition or weights for %s'...
                            ' cubature are not correct.\n'...
                            'Area from section-properties was %.3G'...
                            ' and got %.3G\nError-%% is %.2G'...
                            ' which is higher than allowed %.3G\n'],...
                            cub,A,Ac,error_percent,...
                            ao.max_area_error_percent);
                        continue;
                    end
                end
                x_s=cao.spr.sc(1)-cao.spr.c(1);
                y_s=cao.spr.sc(2)-cao.spr.c(2);
                Io=weights'*(w_at_centers.^2);
                Qo=weights'*w_at_centers;
                Ixo=weights'*(cao.centers(:,1).*w_at_centers);
                Iyo=weights'*(cao.centers(:,2).*w_at_centers);
                Iw=Io-Qo^2/A-y_s*Ixo+x_s*Iyo;      
                elapsed=toc;
                cub_with_card=sprintf("%s(%G)",cub,card);
                if anynan(Iw)
                    fprintf("model=%s, cub=%s failed\n",...
                        model,cub_with_card);
                    continue;
                end
                [~,~,epc]=get_maiw(t,f);
                if epc>ao.max_epc
                    fprintf(['model=%s, cub=%s rejected ' ...
                        'due to epc of %.3G\n'],...
                        model,cub_with_card,epc);
                    continue;
                end
                if ao.latex
                    fprintf("%s%s-%s & %.1f %s & %.1f %s\n", ...
                    "\hspace{1cm}",model,cub_with_card,...
                    Iw*1e12,"\(10^{-12}\)",epc,"\\");
                else
                    fprintf("model=%s-%s, Iw=%.3g, cub took %.3G ms\n", ...
                        model,cub_with_card,Iw,elapsed*1000);
                end
                tf=@()w(cao.centers(:,1),cao.centers(:,2)); %#ok<NASGU>
                es=sprintf("o.%s_%s_%d_walltime=timeit(tf);",...
                    model,cub,card);
                eval(es);                
                if bitand(ao.plot,4) && ci==1
                    s=sprintf("%s for %s",model,fn);
                    figure('Name',s);
                    switch fitMethod
                        case 'fit'
                        plot(f,[t.x t.y],t.w);
                        s=sprintf("Iw=%.3g, rsquare=%.4g",Iw,gof.rsquare);
                        title(s);
                        case 'tpaps'
                        tps_plot(f,list(i),t);
                        s=sprintf("Iw=%.3g",Iw);
                        title(s);
                    end
                    axis equal;
                    ax=gca;
                    dz=(max(t.w)-min(t.w))/min([max(t.x) max(t.y)]);
                    ax.DataAspectRatio=[1 1 dz];
                end
                if (Iw>maxIw)
                    fprintf(['Model %s using %s for file %s rejected,'...
                    ' Iw=%.3g, max=%.3g\n'],...
                    model,cub,fn,...
                    Iw,maxIw);
                    continue;
                end
                es=sprintf("o.%s_%s_Iw=Iw;",model,cub);
                eval(es);
                if ci==1
                    es=sprintf("o.%s_fit=f;",model);
                    eval(es);
                    switch fitMethod
                        case 'fit'
                        es=sprintf("o.%s_gof=gof;",model);
                        eval(es);
                        es=sprintf("o.%s_output=output;",model);
                        eval(es);
                    end
                end
            end
        end
    end
    c{i}=o;
end
end

function [XV,YV]=getVertices(r,n_r,h,b,t)
% get vertices for polyshape making L with rounder corner
% r=inner radius
% which can be seen as quarter of RHS
    if r<t
        error('r=%g but it may not be <t=%g',r,t);
    end
    if r>b
        error('r=%g but it may not be >b=%g',r,b);
    end
    if n_r==1
        error('n_r must be >1 if r is defined');
    end
    if n_r==0
        n_r=4;
    end
    r_in=r;
    r_out=r+t;
    np=6+2*(n_r-1);
    pa=zeros(np,2);
    %points=zeros(4+n_r);
    ri=1;
    % outer bottom left radius
    [n,npa]=draw_radius([r_out, r_out],r_out,pi,n_r);
    pa(ri:ri+n-1,:)=npa;
    ri=ri+n;
    % bottom right corner
    pa(ri,:)=[b,0];
    pa(ri+1,:)=[b,t];
    ri=ri+2;
    % inner bottom left radius
    [n,npa]=draw_radius([t+r_in, t+r_in],r_in,...
        1.5*pi,n_r,false);
    pa(ri:ri+n-1,:)=npa;
    ri=ri+n;
    % line of symmetry
    pa(ri,:)=[t,h];
    pa(ri+1,:)=[0,h];
    XV=pa(:,1);
    YV=pa(:,2);
end