function save_interpolation_results(f,width,height,oms,filename)
%{
Save results in required format 
x0 and y0 are coordinates used in section properties analysis
aw is w from analytical solution
%}
    full_name=sprintf("gen/%s",filename);
    xv=linspace(0,width,oms);
    yv=linspace(0,height,oms);
    [xm,ym]=meshgrid(xv,yv);
    x0=reshape(xm,[],1);
    y0=reshape(ym,[],1);
    w=feval(f,x0,y0);
    tf=isnan(w);
    if any(tf)
        error("nans during function evaluation")
    end
    x=x0-width/2;
    y=y0-height/2;
    awf=@(x,y)rect_psi(x,y,0,0,6,width/2,height/2);
    aw=awf(x,y);
    max_aw=max(abs(aw));
    error_estimate_a=(w+aw)/max_aw;
    error_estimate_b=(w-aw)/max_aw;
    if sum(abs(error_estimate_a))>sum(abs(error_estimate_b))
        error_estimate_source=error_estimate_b;
    else
        error_estimate_source=error_estimate_a;
    end
    error_percent=round(100*error_estimate_source,2);
    max_error_percent=max(abs(error_percent));
    %T=table(x,y,w,x0,y0,aw,error_percent);% for debugging
    T=table(x,y,w);
    writetable(T,full_name);
    fprintf("Saved %s, max error percent=%.2g\n",...
        full_name,max_error_percent);
end