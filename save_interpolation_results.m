function save_interpolation_results(f,width,height,oms,filename)
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
    awf=@(x,y)rect_psi(x,y,0,0,6,height,width);
    aw=awf(x,y);
    T=table(x0,y0,x,y,w,aw);
    writetable(T,full_name);
    fprintf("Saved %s\n",full_name);
end