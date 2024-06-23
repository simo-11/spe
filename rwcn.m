function [v]=rwcn(n,H,W)
% Calculate constant cn
%{
for Eq 33 in
https://www.researchgate.net/publication/361446204_Efficient_modeling_and_order_reduction_of_new_3D_beam_elements_with_warping_via_absolute_nodal_coordinate_formulation
%}
k1=-1^((n+1)/2);
k2=(8*H^2)/(n^3*pi^3*cosh(n*pi*W/(2*H)));
v=k1*k2;
end