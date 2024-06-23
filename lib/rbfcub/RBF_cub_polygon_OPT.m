
function  [weights,cpus,Vcond,moms,res2,phi_str]=RBF_cub_polygon_OPT(P,centers,...
    RBF_type,RBF_scale,strcub)

%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% This routine computes the weights of a cubature rule on a polygonal
% domain described by a polyshape object "P", whose nodes are the RBF
% "centers", RBF defined by "RBF_type" (see "RBF.m") and each RBF centered
% in "centers(i,:)" has RBF scale "RBF_scale(i)".
% Important: the centers must not be on the boundary of the polygonal
%            region.
%--------------------------------------------------------------------------
% Input:
%--------------------------------------------------------------------------
% P: polygonal domain as polyshape object.
%
% centers: RBF centers as a N x 2 matrix.
%
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
% 9: phi=@(r) (max(0,(1-r))).^2;           % Wendland W0
% 10: phi=@(r) (35*r.^2+18*r+3).*(max(0,(1-r))).^6;         % Wendland W4
% 11: phi=@(r) (32*r.^3+25*r.^2+8*r+1).*(max(0,(1-r))).^8;  % Wendland W6
% 12: phi=@(r) (sqrt(2)/(3*sqrt(pi))*(3*r^2*log(r/(1+sqrt(1-r.^2)))+...
%            (2*r^2+1).*sqrt(1-r^2)))*max(0,(1-r));  % Missing Wendland
% 13 phi=@(r) exp(-r);                     % Matern beta_1=3/2.
% 14 phi=@(r) (1+r).*exp(-r);              % Matern beta_2=5/2.
% 15: phi=@(r) (max(0,(1-r)))^10*(429*r^4 + 450*r^3 + 210*r^2 + 50*r + 5)
%                                           % Wendland W8
% Note that Wendland W8 is available only using the triangulation based
% method.
%
% RBF_scale: a positive scalar: all the RBFs have the same scale.
%            a vector of the same number of rows of centers: scale of the
%            RBF for each center.
%
% res2: it is the moment residual in norm 2, i.e.
%                      norm(moms-RBF_cubmat*weights)
%--------------------------------------------------------------------------
% OUTPUT.
%--------------------------------------------------------------------------
% weights: cubature weights (w.r.t. the nodes "centers").
% cpus: some cputimes:
%       cpus(1): computation of "(cubature) Vandermonde matrix";
%       cpus(2): computation of "moments";
%       cpus(3): computation of "polynomial moments";
%       cpus(4): computation of "weights" (by solving a linear system, via
%                "\");
% Vcond: estimate of the conditioning of the (cubature) Vandermonde matrix.
% moms: RBF moments;
% res2: it is "res2=norm(moms-RBF_cubmat*weights,2)" i.e. residuals moments 
%       matching in the 2 norm.
%--------------------------------------------------------------------------%--------------------------------------------------------------------------
% Reference papers:
% 1. R. Cavoretto, A. De Rossi, A. Sommariva, M. Vianello
% RBFCUB: a numerical package for near-optimal meshless cubature on general
% polygons
%
% 2. A. Sommariva and M. Vianello,
% RBF moment computation and meshless cubature on general polygonal regions.

%--------------------------------------------------------------------------
% Copyrights.
%--------------------------------------------------------------------------
%% Copyright (C) 2021- R. Cavoretto, A. De Rossi, A. Sommariva, M. Vianello
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
%% Authors:  
%%          Roberto Cavoretto <roberto.cavoretto@unito.it>
%%          Alessandra De Rossi   <?alessandra.derossi@unito.it>
%%          Alvise Sommariva <alvise@euler.math.unipd.it>
%%          Marco Vianello   <marcov@euler.math.unipd.it>
%%
%% Last Update: July 9, 2021.
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% TROUBLESHOOTING.
%--------------------------------------------------------------------------

warning off;

weights=[]; cpus=[];
if nargin < 1, error('\n \t Polygon not provided'); return; end
if nargin < 2, error('\n \t RBF centers not provided'); return; end
if nargin < 3, RBF_type=5; end
if nargin < 4, RBF_scale=ones(size(centers,1),1); end
if nargin < 5
    if RBF_type <= 3 || RBF_type >= 15
        strcub='triangulation';  
    else 
        strcub='polar'; 
    end
end
if strcmp(strcub,'')
    if  RBF_type <= 3 || RBF_type >= 15
        strcub='triangulation';  
    else
        strcub='polar'; 
    end
end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% MAIN PROGRAM STARTS HERE
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% ...................... building integration matrix ......................

[phi,phi_str] = RBF(RBF_type);

switch RBF_type
    case 1, RBF_order=1;
    case 5, RBF_order=2;
    case 6, RBF_order=2;
    case 7, RBF_order=3;
    case 8, RBF_order=4;
    otherwise, RBF_order=0;
end

% ... adding polynomial components (if required) ...
tic;
Npts=size(centers,1);
switch RBF_order
    case 1
        border_RBF_cubmat=ones(Npts,1);
    case 2
        border_RBF_cubmat=[ones(Npts,1)  centers];
    case 3
        c1=(centers(:,1)); c2=(centers(:,2));
        border_RBF_cubmat=[ones(Npts,1) centers c1.^2 c1.*c2 c2.^2];
    case 4
        c1=(centers(:,1)); c2=(centers(:,2));
        border_RBF_cubmat=[ones(Npts,1) centers c1.^2 c1.*c2 c2.^2 ...
            c1.^3 (c1.^2).*c2 c1.*(c2.^2) c2.^3];
end

% ... computing integration matrix (not the interpolation one) ...
RBF_scales_mat=repmat(1./RBF_scale,1,Npts);
DM = points2distances(centers);
RBF_cubmat=phi(DM.*RBF_scales_mat);

if RBF_order >= 1
    zeromat=zeros(size(border_RBF_cubmat,2));
    RBF_cubmat=[RBF_cubmat border_RBF_cubmat; ...
        border_RBF_cubmat' zeromat];
end

% ......................... computing RBF moments .........................

tic;
if strcmp(strcub,'polar')
    moms=RBF_moms_polar(P,centers,RBF_type,RBF_scale);
else
    moms=RBF_moms_tri(P,centers,RBF_type,RBF_scale);
end
cpus(2)=toc;


% ... adding polynomial terms (if required) ...
tic;
switch RBF_order
    case 1
        [xyw,xvc,yvc,pgon,tri]=polygauss_2018(1,P);
        XX=xyw(:,1); YY=xyw(:,2); WW=xyw(:,3);
        moms(Npts+1,1)=sum(WW);
    case 2
        [xyw,xvc,yvc,pgon,tri]=polygauss_2018(1,P);
        XX=xyw(:,1); YY=xyw(:,2); WW=xyw(:,3);
        moms(Npts+1,1)=sum(WW);
        moms(Npts+2,1)=XX'*WW;
        moms(Npts+3,1)=YY'*WW;
    case 3
        [xyw,xvc,yvc,pgon,tri]=polygauss_2018(2,P);
        XX=xyw(:,1); YY=xyw(:,2); WW=xyw(:,3);
        moms(Npts+1,1)=sum(WW);
        moms(Npts+2,1)=XX'*WW;
        moms(Npts+3,1)=YY'*WW;
        moms(Npts+4,1)=(XX.^2)'*WW;
        moms(Npts+5,1)=(XX.*YY)'*WW;
        moms(Npts+6,1)=(YY.^2)'*WW;
    case 4
        [xyw,xvc,yvc,pgon,tri]=polygauss_2018(3,P);
        XX=xyw(:,1); YY=xyw(:,2); WW=xyw(:,3);
        moms(Npts+1,1)=sum(WW);
        moms(Npts+2,1)=XX'*WW;
        moms(Npts+3,1)=YY'*WW;
        moms(Npts+4,1)=(XX.^2)'*WW;
        moms(Npts+5,1)=(XX.*YY)'*WW;
        moms(Npts+6,1)=(YY.^2)'*WW;
        moms(Npts+7,1)=(XX.^3)'*WW;
        moms(Npts+8,1)=(XX.^2.*YY)'*WW;
        moms(Npts+9,1)=(XX.*YY.^2)'*WW;
        moms(Npts+10,1)=(YY.^3)'*WW;
end
cpus(3)=toc;

% ......................... computing RBF weights .........................

ls_meth=2;

tic;
switch ls_meth
        case 0 % BETTER SOLUTION.
        % polynomial basis orthogonalization
        V=RBF_cubmat'; m=moms;
        [Q,R]=qr(V,0); orthmom=(m'/R)'; 
        weights0=(V/R)'\orthmom;
    case 1
        % polynomial basis orthogonalization
        [Q,R]=qr(RBF_cubmat',0); Q=real(Q); orthmom=R'\moms; 
        weights0=Q'\orthmom;
    otherwise
        % no polynomial basis orthogonalization
        weights0=RBF_cubmat\moms;
end
cpus(1,4)=toc;

weights=weights0;
switch RBF_order
    case 1, weights=weights(1:end-1);
    case 2, weights=weights(1:end-3);
    case 3, weights=weights(1:end-6);
    case 4, weights=weights(1:end-10);
end

if nargout >= 3, Vcond=cond(RBF_cubmat); end
if nargout >= 5, res2=norm(moms-RBF_cubmat*weights0,2);  end

% +++++++++++++++++++++++++ RBF_cub ENDS HERE +++++++++++++++++++++++++++++







%==========================================================================
% ATTACHED FUNCTIONS.
%==========================================================================

function distances = points2distances(points)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% Create distance matrix "distances" w.r.t. "points", i.e.
%      distances(i,j)=distance(points(i,:)-points(j,:)).
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% points: set of points whose coordinates are stored as N x 2 matrix.
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% distances: distance matrix, where
%      distances(i,j)=distance(points(i,:)-points(j,:)).
%--------------------------------------------------------------------------

% Get dimensions.
[card,dim]=size(points);

% All inner products between points.
distances=points*points';

% Vector of squares of norms of points.
lsq=diag(distances);

% Distance matrix.
distances=sqrt(repmat(lsq,1,card)+repmat(lsq,1,card)'-2*distances);















