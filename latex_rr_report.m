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
n=11;
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
fig=figure(317);
fp=sprintf("Selected points for shear stresses at corner");
fig.Name=fp;
plot(ps);
axis equal;
for vi=1:size(xvalues,2)
    x_s=xvalues(vi)-cao.spr.sc(1);
    y_s=yvalues(vi)-cao.spr.sc(2);
    txt=sprintf("(%.3G,%.3G)",x_s,y_s);
    text(xvalues(vi),yvalues(vi),txt);
end
for mi=1:ms
    model=ao.models(mi);
    fprintf(tlFileID,"%s%s\n",...
        model,"\\");
    es=sprintf("f=ro.%s_fit;",model);
    eval(es);
    w=@(x,y)f(x,y);
    wvalues=w(xvalues,yvalues);
    [dxvalues,dyvalues]=differentiate(f,xvalues,yvalues);
    for vi=1:size(xvalues,2)
        x_s=xvalues(vi)-cao.spr.sc(1);
        y_s=yvalues(vi)-cao.spr.sc(2);
        fprintf(tlFileID, ...
            "%.3G & %.3G & %.3G & %.3G & %.3G%s\n",...
        x_s,y_s,wvalues(vi),...
        dxvalues(vi),dyvalues(vi),"\\");
    end
end
fprintf("Wrote %s\n",tlfn);
fclose(tlFileID);    
