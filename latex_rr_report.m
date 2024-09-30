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
fl=1.6; % fraction of length added after curve >sqrt(2)
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
[t,~,~]=define_scattered_pointset(30,domain);
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
        qr=quiver(xvalues,yvalues,-yvalues,xvalues,0.8);
        qr.Color="red";
        scale=qr.ScaleFactor;
        qw=quiver(xvalues,yvalues,scale*dxvalues,scale*dyvalues,'off');
        qw.Color="blue";        
    end
    fig=figure(317+mi);
    hold off;
    fp=sprintf("Shear stresses at corner using %s",model);
    fig.Name=fp;
    plot(ps);
    axis equal;
    hold on;
    [dxvalues,dyvalues]=differentiate(f,X,Y);
    % rotation related values first as they are larger
    qr=quiver(X,Y,-Y,X,0.8);
    qr.Color="red";
    scale=qr.ScaleFactor;
    qw=quiver(X,Y,scale*dxvalues,scale*dyvalues,'off');
    qw.Color="blue";
    fn=sprintf("rr_shear_stress_using_%s",model);
    save_pdf(fig,fn);
end
fprintf("Wrote %s\n",tlfn);
fclose(tlFileID);
save_pdf(pfig,"points_for_stress");
