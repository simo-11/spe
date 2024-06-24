function [maiw,im,epc]=get_maiw(t,f)
% gets maximum absolute value of interpolated function w mainw
% im abs of maximum input 
% epc error per cent
    [M,I]=max(t.w);
    [m,i]=min(t.w);
    if(abs(M)>abs(m))
        x=t.x(I);
        y=t.y(I);
        im=t.w(I);
    else
        x=t.x(i);
        y=t.y(i);
        im=t.w(i);
    end
    if isa(f,'sfit')
        maiw=f(x,y);
    else 
        maiw=fnval(f,[x y]');
    end
    im=abs(im);
    maiw=abs(maiw);
    epc=100*(abs(maiw-im))/im;
end