%% settings
ao.models=["poly44","cubicinterp","tps"];
ao.cubs=["integral2","glaubitz","rbfcub"];
ao.u_models=["poly44","cubicinterp","tps"];
ao.u_cubs=["rbfcub"]; %#ok<NBRAK2>
ao.rhs_models=["poly44","cubicinterp","thinplateinterp"];
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
E=210e9;
G=E/2.6;
L=0.5;
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
%% SHS 150x150x8 with rounded corners using n_r=8
rr8=testRHS(height=150,width=150,t=8,r=16,n_r=8,models=ao.rhs_models,...
    cubs=ao.rhs_cubs,debug=ao.debugLevel,cards=ao.rhs_cards);
%% latex report for first SHS 150x150x8 with rounded corners, n_r=8
% > echodemo('paper_cells',10)
latex_rr_report(rr8{1})
%% SHS 150x150x8 with rounded corners using n_r=24
rr24=testRHS(height=150,width=150,t=8,r=16,n_r=24,models=ao.rhs_models,...
    cubs=ao.rhs_cubs,debug=ao.debugLevel,cards=ao.rhs_cards);
%% latex report for first SHS 150x150x8 with rounded corners, n_r=24
% > echodemo('paper_cells',12)
rr24{1}.ao.rr_card=150; % 
rr24{1}.ao.rr_fl=0; % ~1.6 for paper, max 12.3
latex_rr_report(rr24{1})
%% differential equation for torsion 
syms theta(z) git eiw T
ode=git*diff(theta,z)-eiw*diff(theta,z,3)==T;
Dz=diff(theta);
Dz2=diff(theta,2);
c1=theta(0)==0;
c2=Dz(0)==0;
c3=Dz2(L)==0;
conds=[c1 c2 c3];
sol(z)=dsolve(ode,conds);
[ssol,k]=subexpr(sol,k);
%% Plot distribution of rotation and derivates in z-direction
fp=sprintf("gen/results*.json");
list=dir(fp);
n=size(list,1);
prf_count=5;
c=cell(prf_count,1);
ci=0;
for i=1:n
    name=list(i).name;
    prf=regexp(name,'results-(?<name>.*)-(?<nodes>\d+).json','names');
    fn=sprintf("%s/%s",list(i).folder,name);
    spr=jsondecode(fileread(fn));
    % pick good ones
    if contains(name,"box-175-129")
        key="rhs";
    elseif contains(name,"cold-formed-u-100-50-4-8-12")
        key="u";
    elseif contains(name,"rectangle-100-100-1317")
        key="r100";
    elseif contains(name,"rectangle-30-100-347")
        key="r30";
    elseif contains(name,"rhs-150-150-8-16-8-1716")
        key="shs";
    else
        continue
    end
    o.key=key;
    o.spr=spr;
    o.prf=prf;
    o.k=sqrt(G*spr.j/(E*spr.gamma));
    fprintf("%s & %.3G & %.3G& %.3G%s\n",...
        prf.name,spr.gamma*1e12,spr.j*1e6,o.k,"\\");
    ci=ci+1;
    c{ci}=o;
end
x = linspace(0,L);
fig=figure(313);
theta_plot(0,c,x,prf_count,L);
fig.Name="rotation";
save_pdf(fig,fig.Name);
fig=figure(fig.Number+1);
theta_plot(1,c,x,prf_count,L);
fig.Name="warping";
save_pdf(fig,fig.Name);
fig=figure(fig.Number+1);
theta_plot(2,c,x,prf_count,L);
fig.Name="bimoment";
save_pdf(fig,fig.Name);
function theta_plot(d_level,c,x,prf_count,L)
    clf;
    hold on;
    for i=1:prf_count
        k=c{i}.k;
        key=c{i}.key;
        switch d_level
            case 0
                y0=tanh(k*L)*(cosh(k*x)-1)-sinh(k*x)+k*x;
                ytext='rotation/max rotation';
            case 1
                y0=tanh(k*L)*sinh(k*x)-cosh(k*x)+1;
                ytext='warping/max warping';
            case 2
                y0=tanh(k*L)*cosh(k*x)-sinh(k*x);
                ytext='bimoment/max bimoment';
            otherwise
                error("d_level %d is not pupported",d_level)
        end
        max_y=max(y0);
        y=y0./max_y;
        plot(x,y,'DisplayName',key);
    end
    legend(Location="southeast")
    xlabel('z[m]');
    ylabel(ytext);
end

