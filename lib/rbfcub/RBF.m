function [phi,phi_str] =  RBF (RBF_type)

%--------------------------------------------------------------------------
% Input:
%--------------------------------------------------------------------------
% RBF_type: the choosen RBF is the "RBF_type"-th in the list of functions
%     described in "RBF.m".
%
% 1: phi=@(r) (1+r.*r).^(1/2);             % Multiquadric
% 2: phi=@(r) exp(-r.*r);                  % Gaussian
% 3: phi=@(r) (1+r.*r).^(-1/2);            % Inverse Multiquadric
% 4: phi=@(r) (1+4*r).*(max(0,(1-r))).^4;  % Wendland 2
% 5: phi=@(r) r.*r.*log(r+eps*(1-sign(r))); % TPS
% 6: phi=@(r) r.^3;                        % polyharmonic spline
% 7: phi=@(r) r.^5;                        % polyharmonic spline
% 8: phi=@(r) r.^7;                        % polyharmonic spline
% 9: phi=@(r) (max(0,(1-r))).^2;             % Wendland W0
% 10: phi=@(r) (35*r.^2+18*r+3).*(max(0,(1-r))).^6;         % Wendland W4
% 11: phi=@(r) (32*r.^3+25*r.^2+8*r+1).*(max(0,(1-r))).^8;  % Wendland W6
% 12: phi=@(r) (sqrt(2)/(3*sqrt(pi))*(3*r^2*log(r/(1+sqrt(1-r.^2)))+...
%            (2*r^2+1).*sqrt(1-r^2)))*max(0,(1-r));  % Missing Wendland
% 13 phi=@(r) exp(-r);                     % Matern beta_1=3/2.
% 14 phi=@(r) (1+r).*exp(-r);              % Matern beta_2=5/2.
% 15: phi=@(r) (max(0,(1-r)))^10*(429*r^4 + 450*r^3 + 210*r^2 + 50*r + 5)
%--------------------------------------------------------------------------
% Output:
%--------------------------------------------------------------------------
% phi: RBF function.
% phi_str: string with the RBF "name".
%--------------------------------------------------------------------------
% Copyrights.
%--------------------------------------------------------------------------
%% Copyright (C) 2007-2019 Alvise Sommariva, Marco Vianello.
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%%
%% Author:  Alvise Sommariva <alvise@euler.math.unipd.it>
%%          Marco Vianello   <marcov@euler.math.unipd.it>
%%
%% Date: February 14, 2020.
%--------------------------------------------------------------------------

switch RBF_type
    case 1
        phi_str='MQ';
        phi=@(r) (1+r.*r).^(1/2);            % Multiquadric
    case 2
        phi_str='Gaussian';
        phi=@(r) exp(-r.*r);                  % Gaussian
    case 3
        phi_str='IMQ';
        phi=@(r) (1+r.*r).^(-1/2);            % Inverse Multiquadric
    case 4
        phi_str='W2';
        phi=@(r) (1+4*r).*(max(0,(1-r))).^4;  % Wendland 2
    case 5
        phi_str='TPS';
        pert=@(r) eps*(1-abs(sign(r)));
        phi=@(r) r.*r.*log(r+pert(r)); % TPS
    case 6
        phi_str='PH: r^3';
        phi=@(r) r.^3;                        % polyharmonic spline
    case 7
        phi_str='PH: r^5';
        phi=@(r) r.^5;                        % polyharmonic spline
    case 8
        phi_str='PH: r^7';
        phi=@(r) r.^7;                        % polyharmonic spline
    case 9
        phi_str='W0';
        phi=@(r) (max(0,(1-r))).^2;             % Wendland W0
    case 10
        phi_str='W4';
        phi=@(r) (35*r.^2+18*r+3).*(max(0,(1-r))).^6;  % Wendland W4
    case 11
        phi_str='W6';
        phi=@(r) (32*r.^3+25*r.^2+8*r+1).*(max(0,(1-r))).^8;  % Wendland W6
    case 12                                    % Missing Wendland
        phi_str='Missing Wendland';
        pert=@(r) eps*(1-abs(sign(r)));
        phi=@(r) (sqrt(2)/(3*sqrt(pi))*...
            (3*r.^2.*log( r./(1+sqrt(1-r.^2)) + pert(r) )+...
            (2*r.^2+1).*sqrt(1-r.^2))).*max(0,(1-r));
    case 13                             % Matern beta_1=(d+1)/2, where d=2.
        phi_str='Matern beta_1';
        phi=@(r) exp(-r);
    case 14                             % Matern beta_2=(d+3)/2, where d=2.
        phi_str='Matern beta_2';
        phi=@(r) (1+r).*exp(-r);
    case 15
        phi_str='Wendland W8';
        phi=@(r) (max(0,(1-r))).^10.*(429*r.^4 + 450*r.^3 + 210*r.^2 + ...
            50*r + 5);
        
    otherwise
        error('RBF type not implemented')
end




