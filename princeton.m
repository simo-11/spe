%% settings
%ao.models=["poly44","cubicinterp","tps"];
ao.models=["cubicinterp"]; %#ok<*NBRAK2>
%ao.cubs=["integral2","glaubitz","rbfcub"];
ao.cubs=["integral2"];
ao.debugLevel=0;
ao.plot=0;
ao.scat_type='halton';
ao.cards=[50,100,400,500];
ao.rsquareMin=0.9;
ao.n="294";%number of nodes in section properties model
%ao.n="375";
ao.oms=10;% output mesh size, number of points in each direction
add_lib_to_path
E=71.7e9;
nu=0.31;
G=E/(2*(1+nu));
L=0.508;
w=3.2024;% mm
h=12.377;% mm
%w=100;
%h=10;
width=w/1000;% m
height=h/1000;% m
%% analytical solution
thin_value=width^3*height^3/144;
fprintf("Analytical Iw for %Gx%G, for thin: %.3G\n",w,h, thin_value);
for i=6:-1:0
    tic;
    n=floor(i);
    Iw=iw_rect(height=height,width=width,n=n);
    elapsed=toc;
    fprintf("Iw=%.3G (n=%G) took %G ms\n",Iw,n,elapsed*1000);
end
%% interpolate based on section properties results
expected_name=sprintf("gen/warping-rectangle-%g-%g-%s.csv",h,w,ao.n);
if ~isfile(expected_name)
    error('results are not available, expected to have file %s',...
        expected_name);
end
r=testRectangle(height=h,width=w,models=ao.models,...
    cubs=ao.cubs,debug=ao.debugLevel,cards=ao.cards,n=ao.n,...
    latex=0);
if size(r,1)>0
    r{1} %#ok<NOPTS>
end
    error('')
%% write warping results based on interpolation f
list=dir(r{1}.file);
expected_name=sprintf("warping-rectangle-%g-%g-%s.csv",h,w,ao.n);
if list.name~=expected_name
    error('results are for %s,\nbut current model expects %s',...
        list.name,expected_name);
end
ms=size(ao.models,2);
ao.oms=4;
for mi=1:ms
    model=ao.models(mi);
    es=sprintf("f=r{1}.%s_fit;",model);
    eval(es);
    es=sprintf("fn='%s-%g-%g-%s-%dx%d.xlsx';",...
        model,h,w,ao.n,ao.oms,ao.oms);
    eval(es);
    save_interpolation_results(f,width,height,ao.oms,fn);
end
