
function [w,dmax]=glaubitz_algorithm(xS,domain_structure,wQMC,approach,...
    d_start)

%--------------------------------------------------------------------------
% INPUT
%--------------------------------------------------------------------------
% xS: scattered data
% domain_structure: structure of the domain; it must include:
%       domain.ratio=area_domain/area_bounding_box;
% wQMC: weights of a QMC cubature rule, having "t" as nodes;
% approach: 'LS' or 'l1'; 
% d_start : starting degree of analysis;
%
% The variables "approach" and "d_start" are not mandatory.
%--------------------------------------------------------------------------
% OUTPUT
%--------------------------------------------------------------------------
% w : weights of the cubature rule, having "t" as nodes;
% d : highest possible degree of exactness (ADE of (t,w))
%--------------------------------------------------------------------------
% SUBROUTINE USED
%--------------------------------------------------------------------------
% 1. compute_weights (attached)
% 2. define monomials (attached)
% 3. define_cub_rule (external, compute an alg. rule, with prescribed 
%    degree of exactness, in several domains)
%--------------------------------------------------------------------------

% ................... troubleshooting ................... 

if nargin < 4, approach='LS'; end
if isempty(approach), approach='LS'; end
if not(strcmp(approach,'LS') | strcmp(approach,'l1')), approach='LS'; end

if nargin < 5, d_start=0; end
if isempty(d_start), d_start=0; end


% ................ assign some variables ................ 

xS_N=size(xS,1);
xS_dim=size(xS,2);
R = diag(wQMC); % discrete weights matrix
w = zeros(xS_N,1); % cubature weights


% ................ set up initial values ................ 

d = d_start; % degree of exactness (DoE)
K = nchoosek(xS_dim + d, xS_dim); % number of basis elements
check_rank = K; % rank
w_min = min(w); % smallest cubature weight

w_history={};

% ....... loop in which d is increased until the CF becomes unstable ......

e = -1e-14; % tollerance to allow small rounding errors
while check_rank == K && w_min >= e

    d = d+1; % increase the DoE
    K = nchoosek(xS_dim + d, xS_dim); % number of basis elements

    basis=define_monomials(d,xS_dim);
    P = basis(xS); % Vandermonde matrix

    check_rank = rank(P); % rank of P

    % compute moments
    [XW,~,domain_structure]=define_cub_rule(domain_structure,d);

    PXW=basis(XW(:,1:2));
    m=PXW*XW(:,3);

    % check wheter the data points are d-unisolvent
    if check_rank == K

        % different approaches
        switch approach
            case 'LS'
                A = P*sqrt(R);
                v = lsqminnorm(A,m); % indirect computation by optim. tools
                w = sqrt(R)*v;
            case 'l1'
                options = optimoptions('linprog','Display','none');
                w = linprog(ones(xS_N,1), [], [], P, m, zeros(xS_N,1),...
                    [], options ); % l1 weights
        end

    end

    if sum(w)==sum(w) && sum(abs(w))~=Inf && ...
            length(w) > 0 % w does not contain NaNs or Infs and is nonempty
        w_min = min(w);
    else
        w_min = -1;
    end

    w_history{end+1}=w;

end

w=w_history{end-1};
dmax=d-1;






function basis=define_monomials(d,dim)

K = nchoosek(dim + d, dim); % binomial coefficient/ dimension
alpha1 = zeros(K,1); % vector of exponents for x
alpha2 = zeros(K,1); % vector of exponents for y
alpha3 = zeros(K,1); % vector of exponents for z
m = zeros(K,1); % moments

%% exponents and basis
if dim == 1
    alpha1 = (0:d)';
    basis = @(x) x'.^alpha1; % basis
elseif dim == 2
    k = 1;
    for k1=0:1:d
        for k2=0:1:d
            if k1+k2<=d
                alpha1(k) = k1;
                alpha2(k) = k2;
                k = k+1;
            end
        end
    end
    basis = @(x) x(:,1)'.^alpha1 .* x(:,2)'.^alpha2; % basis
elseif dim == 3
    k = 1;
    for k1=0:d
        for k2=0:d
            for k3=0:d
            	if k1+k2+k3<=d
                	alpha1(k) = k1;
                    alpha2(k) = k2;
                    alpha3(k) = k3;
                    k = k+1;
                end
            end
        end
    end
    basis = @(x) x(:,1)'.^alpha1 .* x(:,2)'.^alpha2 .* x(:,3)'.^alpha3; % basis
else
    error('Desired dimension not yet implemented!')
end










