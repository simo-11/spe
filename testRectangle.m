function [c]=testRectangle(ao)
arguments
    ao.height=100
    ao.width=100
    ao.models=["poly44","cubicinterp","tps"]
    ao.cubs=["integral2","glaubitz","rbfcub"]
    ao.scat_type='halton'
    ao.cards=[500]
    ao.rsquareMin=0.9
    ao.debugLevel=0
    ao.plot=0
    ao.latex=1
    ao.n="*"
end
%{
Test warping function fit using csv files in gen directory
%}
fp=sprintf("gen/warping-rectangle-%g-%g-%s.csv",...
    ao.height,ao.width,ao.n);
list=dir(fp);
n=size(list,1);
c=cell(n,1);
ms=size(ao.models,2);
cs=size(ao.cubs,2);
H=ao.height/1000;
W=ao.width/1000;
domain.vertices=[0 0; W 0; W H; 0 H; 0 0];
XV=domain.vertices(:,1);
YV=domain.vertices(:,2);
domain.polyshape=polyshape(XV,YV);
domain.domain='rectangle';
[xlimit,ylimit]=boundingbox(domain.polyshape);
domain.dbox=[xlimit; ylimit];
cao.debugLevel=ao.debugLevel;
cao.plot=ao.plot;
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
                w=@(x,y)f(x,y).^2;
            case 'tpaps'
                if strlength(model)>3
                    pin=str2double(extractAfter(model,3));
                else
                    pin=1;
                end
                f=tpaps([t.x t.y]',t.w',pin);
                w=@(x,y)reshape(fnval(f,[x(:)';y(:)']).^2,...
                    size(x,1),[]);
        end
        if ao.debugLevel>1
            disp(f);
            switch fitMethod
                case 'fit'
                disp(gof);
            end
        end
        for ci=1:cs
            cub=ao.cubs(ci);
            cub_with_card="";
            for card=ao.cards
                if cub==cub_with_card
                    continue
                end
                cao.card=card;
                [cao.centers,cao.dbox,cao.area_domain]=...
                   define_scattered_pointset(card,domain,ao.scat_type);            
                cao.w_at_centers=w(cao.centers(:,1),cao.centers(:,2));
                tic
                [Iw,cub_suffix]=do_cub(w,domain,cub,cao);
                elapsed=toc;
                cub_with_card=sprintf("%s%s",cub,cub_suffix);
                if anynan(Iw)
                    fprintf("model=%s, cub=%s failed\n",...
                        model,cub_with_card);
                    continue;
                end
                if ao.latex
                    fprintf("%s%s-%s & %.3g %s\n", ...
                    "\hspace{1cm}",model,cub_with_card,...
                    Iw*1e12,"\(10^{-12}\)\\");
                else
                    fprintf("model=%s-%s, Iw=%.3g, cub took %.3G ms\n", ...
                    model,cub_with_card,Iw,elapsed*1000);
                end
                if ao.plot && ci==1
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


