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
ds2=(R+T)/sqrt(2);
xvalues=[W-R W-R W-ds W-ds2 W   W-T];
yvalues=[H   H-T H-ds H-ds2 H-R H-R];
%{
pull points a little bit out of bounds to allow
numerical differentiation to be done without extrapolation
as is needed in cubicinterp
%}
mx = 2*eps^(1/3)*max(1,W);
my = 2*eps^(1/3)*max(1,H);
xvalues=xvalues-mx;
yvalues=yvalues-my;
rfn=replace(ro.file,"warping","results");
rfn=replace(rfn,"csv","json");
cao.spr=jsondecode(fileread(rfn));
tlfn=replace(ro.file,"warping","latex-report");
tlfn=replace(tlfn,".csv",".ltx");
tlFileID=fopen(tlfn,'w');
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
