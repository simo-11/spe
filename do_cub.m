function [I,cub_suffix]=do_cub(f,domain,cub,ao)
%
%{
f=function to evaluate
domain=struct with fields
 polyshape
 dbox
 vertices
cub=cubature method
t=table with columns from section-properties
 x,y
 w warping
ao=argument object struct with fields
 debugLevel int
 centers matrix of size ncenters,2
 w_at_centers
%}
arguments
    f
    domain
    cub
    ao
end
if ao.plot
    plot(domain.polyshape);
    hold on;
    axis equal;
    plot(ao.centers(:,1),ao.centers(:,2),'o')
    hold off
end
switch cub
    case 'integral2'
        I=integral2(f,domain.dbox(1,1),domain.dbox(1,2),...
            domain.dbox(2,1),domain.dbox(2,2));
        if nargout>1
           cub_suffix="";
        end
        return
    otherwise
        w=get_weights(domain,cub,ao);
end
if nargout>1
   cub_suffix=sprintf("(%G)",ao.card); 
end
I=w'*ao.w_at_centers;
end