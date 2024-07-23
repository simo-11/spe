%% settings
ao.models=["poly44","cubicinterp","tps"];
ao.cubs=["integral2","glaubitz","rbfcub"];
ao.u_models=["poly44","cubicinterp","tps"];
ao.u_cubs=["rbfcub"]; %#ok<NBRAK2>
ao.rhs_models=["poly44","cubicinterp","tps"];
ao.rhs_cubs=["rbfcub"]; %#ok<NBRAK2>
ao.debugLevel=0;
ao.plot=0;
ao.scat_type='halton';
ao.cards=[50,100,400,500];
ao.u_cards=[50,100,200,400];
ao.rhs_cards=[50,100,200,400];
ao.rsquareMin=0.9;
ao.n="*";
add_lib_to_path
%% analytical square solid rectangle 100x100
w=100;
h=100;
width=w/1000;
height=h/1000;
fprintf("Analytical Iw for %Gx%G\n",w,h);
for i=0:6
    tic;
    Iw=iw_rect(height=height,width=width,n=i);
    elapsed=toc;
    fprintf("Iw=%.3G (n=%G) took %G ms\n",Iw,i,elapsed*1000);
end
%% analytical solid rectangle 100x10
w=100;
h=10;
fprintf("Analytical Iw for %Gx%G\n",w,h);
width=w/1000;
height=h/1000;
thin_value=width^3*height^3/144;
fprintf("Analytical Iw for %Gx%G, for thin %.3G\n",w,h, thin_value);
for i=6:-1:0
    tic;
    n=floor(i);
    Iw=iw_rect(height=height,width=width,n=n);
    elapsed=toc;
    fprintf("Iw=%.3G (n=%G) took %G ms\n",Iw,n,elapsed*1000);
end
%% square solid rectangle 100x100
ss=testRectangle(height=100,width=100,models=ao.models,...
    cubs=ao.cubs,debug=ao.debugLevel,cards=ao.cards);
%% solid rectangle 100x10
r=testRectangle(height=10,width=100,models=ao.models,...
    cubs=ao.cubs,debug=ao.debugLevel,cards=ao.cards);
%% sharp cornered U-section 100x50x4
su=testU(height=100,width=50,t=4,r=0,n_r=0,models=ao.u_models,...
    cubs=ao.u_cubs,debug=ao.debugLevel,cards=ao.u_cards ...
    );
%% U-section 100x50x4 with rounded corners
ru=testU(height=100,width=50,t=4,r=8,n_r=8,models=ao.u_models,...
    cubs=ao.u_cubs,debug=ao.debugLevel,cards=ao.u_cards);
%% sharp cornered SHS 150x150x8
sr=testRHS(height=150,width=150,t=8,r=0,n_r=0,models=ao.rhs_models,...
    cubs=ao.rhs_cubs,debug=ao.debugLevel,cards=ao.rhs_cards ...
    );
%% SHS 150x150x8 with rounded corners
rr=testRHS(height=150,width=150,t=8,r=16,n_r=8,models=ao.rhs_models,...
    cubs=ao.rhs_cubs,debug=ao.debugLevel,cards=ao.rhs_cards);
