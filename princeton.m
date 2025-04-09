%% settings
%ao.models=["poly44","cubicinterp","tps"];
ao.models=["cubicinterp"];
%ao.cubs=["integral2","glaubitz","rbfcub"];
ao.cubs=["integral2"];
ao.debugLevel=0;
ao.plot=0;
ao.scat_type='halton';
ao.cards=[50,100,400,500];
ao.rsquareMin=0.9;
ao.n="294";
add_lib_to_path
E=71.7e9;
nu=0.31;
G=E/(2*(1+nu));
L=0.508;
w=3.2024;
h=12.377;
%% analytical solid rectangle
width=w/1000;
height=h/1000;
thin_value=width^3*height^3/144;
fprintf("Analytical Iw for %Gx%G, for thin %.3G\n",w,h, thin_value);
for i=6:-1:0
    tic;
    n=floor(i);
    Iw=iw_rect(height=width,width=height,n=n);
    elapsed=toc;
    fprintf("Iw=%.3G (n=%G) took %G ms\n",Iw,n,elapsed*1000);
end
%% solid rectangle
r=testRectangle(height=h,width=w,models=ao.models,...
    cubs=ao.cubs,debug=ao.debugLevel,cards=ao.cards,n=ao.n,...
    latex=0);
%% write warping results based on interpolation f
f=r{1}.cubicinterp_fit;