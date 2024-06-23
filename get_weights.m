function w=get_weights(domain,cub,ao)
%
%{
domain=struct with fields
 polyshape
 dbox
 vertices
cub=cubature method
ao=argument object struct with fields
 debugLevel int
 centers matrix of size ncenters,2
%}
arguments
    domain
    cub
    ao
end
switch cub
    case 'glaubitz'
        area_dbox=diff(ao.dbox(1,:))*diff(ao.dbox(2,:));
        wQMC=(ao.area_domain/area_dbox)/size(ao.centers,1);
        [w,deg_rule]=glaubitz_algorithm(ao.centers,...
            domain,wQMC,'LS');%#ok<ASGLU>
    case 'rbfcub'
        [w, cpus, Vcond , moms , res2 ,phi_str]=...
            RBF_cub_polygon_OPT(domain.polyshape,ao.centers);%#ok<ASGLU>
    otherwise
        error("cub value %s is not supported\n",cub)
end
