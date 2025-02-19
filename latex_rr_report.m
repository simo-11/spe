function latex_rr_report(ro)
%{
produce report for RSH with rounded corners based on result object ro
%}
ao=ro.ao;
ms=size(ao.models,2);
H=ao.height/1000;
W=ao.width/1000;
T=ao.t/1000;
R=ao.r/1000;
ds=R/sqrt(2);
%
% p1  p2
% p6
%          p5
%          p4    p3
% fraction of length added after curve >sqrt(2) and <W-ds
fl=ao.rr_fl;
if fl<=0
    syms nice_fl;
    D=min(W,H);
    nice_fl_value=solve(D-nice_fl*ds==D/2,nice_fl);
    fl=double(nice_fl_value);
end
if fl<sqrt(2)
    error("ro.ao.rr_fl must be >%.3G but was %.3G",sqrt(2),fl);
end
syms max_fl;
max_fl_value=solve(min(W,H)-max_fl*ds==ds,max_fl);
if fl>max_fl_value
    error("ro.ao.rr_fl must be <%.3G but was %.3G",max_fl_value,fl);
end
p1=[W-fl*ds H];
p2=[W-ds H];
p3=[W   H-fl*ds];
p4=[W-T H-fl*ds];
p5=[W-T H-ds];
p6=[W-fl*ds H-T];
cc=[W-R H-R];
n=21;
ci=floor(n/2)+1; % center index
[~,c1]=draw_radius(cc,R,pi/2,n,0,pi/2);
[~,c2]=draw_radius(cc,R-T,0,n,1,pi/2);
points=[p1;p2;c1;p3;p4;p5;c2;p6];
ps=polyshape(points);
%{
values are in upper right corner
pull points abit out of bounds to allow
numerical differentiation to be done without extrapolation
as is needed in cubicinterp
%}
abit=W/1000;
xvalues=[W-R W-R c1(ci,1)-abit c2(ci,1) W-abit   W-T];
yvalues=[H-abit  H-T c1(ci,1)-abit c2(ci,1) H-R H-R];
rfn=replace(ro.file,"warping","results");
rfn=replace(rfn,"csv","json");
cao.spr=jsondecode(fileread(rfn));
tlfn=replace(ro.file,"warping","latex-report");
tlfn=replace(tlfn,".csv",".ltx");
tlFileID=fopen(tlfn,'w');
pfig=figure(317);
hold off;
fp=sprintf("Selected points for shear stresses at corner");
pfig.Name=fp;
plot(ps);
axis equal;
for vi=1:size(xvalues,2)
    x_s=1000*(xvalues(vi)-cao.spr.sc(1));
    y_s=1000*(yvalues(vi)-cao.spr.sc(2));
    txt=sprintf("  (%.3G,%.3G)",x_s,y_s);
    text(xvalues(vi),yvalues(vi),txt);
end
line(xvalues,yvalues,'LineStyle','none','Marker','o')
domain.vertices=ps.Vertices;
domain.polyshape=ps;
domain.domain='rectangle';
[xlimit,ylimit]=boundingbox(domain.polyshape);
domain.dbox=[xlimit; ylimit];
[t,~,~]=define_scattered_pointset(ao.rr_card,domain);
X=t(:,1);
Y=t(:,2);
for mi=1:ms
    model=ao.models(mi);
    fprintf(tlFileID,"%s%s\n",...
        model,"\\");
    es=sprintf("f=ro.%s_fit;",model);
    eval(es);
    w=@(x,y)f(x,y);
    wvalues=w(xvalues,yvalues);
    [dxvalues,dyvalues,fxx,~,fyy]=differentiate(f,xvalues,yvalues);
    for vi=1:size(xvalues,2)
        x_s=1000*(xvalues(vi)-cao.spr.sc(1));
        y_s=1000*(yvalues(vi)-cao.spr.sc(2));
        dx_s=1000*dxvalues(vi);
        dy_s=1000*dyvalues(vi);
        laplace=fxx(vi)+fyy(vi);
        fprintf(tlFileID, ...
            "%.3G & %.3G & %.2f & %.3G & %.3G & %.3G%s\n",...
        x_s,y_s,1e6*wvalues(vi),...
        dx_s,dy_s,laplace,"\\");
    end
    if model=="cubicinterp"
        figure(pfig);
        hold on
        qr=quiver(xvalues,yvalues,...
            cao.spr.sc(2)-yvalues,xvalues-cao.spr.sc(1),0.8);
        qr.Color="red";
        scale=qr.ScaleFactor;
        qw=quiver(xvalues,yvalues,scale*dxvalues,scale*dyvalues,'off');
        qw.Color="blue";
        xlabel("x[m]");
        ylabel("y[m]");
    end
    fig=figure(317+mi);
    hold off;
    fp=sprintf("Shear stresses components at corner using %s",model);
    fig.Name=fp;
    plot(ps);
    axis equal;
    hold on;
    [dxvalues,dyvalues]=differentiate(f,X,Y);
    % rotation related values first as they are larger
    qr=quiver(X,Y,cao.spr.sc(2)-Y,X-cao.spr.sc(1),0.8);
    qr.Color="red";
    scale=qr.ScaleFactor;
    qw=quiver(X,Y,scale*dxvalues,scale*dyvalues,'off');
    qw.Color="blue";
    fn=sprintf("rr_shear_stress_components_using_%s",model);
    xlabel("x[m]");
    ylabel("y[m]");
    save_pdf_and_fig(fig,fn);
    fig=figure(327+mi);
    hold off;
    fp=sprintf("Shear stresses at corner using %s",model);
    fig.Name=fp;
    plot(ps);
    axis equal;
    hold on;
    % rotation related values first as they are larger
    qr=quiver(X,Y,cao.spr.sc(2)-Y+dxvalues,...
        X-cao.spr.sc(1)+dyvalues,0.8);
    qr.Color="red";
    xlabel("x[m]");
    ylabel("y[m]");
    fn=sprintf("rr_shear_stress_using_%s",model);
    save_pdf_and_fig(fig,fn);
end
fprintf("Wrote %s\n",tlfn);
fclose(tlFileID);
save_pdf_and_fig(pfig,"points_for_stress");
