
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