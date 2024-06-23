 
function moms=RBF_moms_polar(P,centers,RBF_type,RBF_scale)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% Once it has been defined the polygonal domain P (in polyshape form),
% this routine computes the moments of a RBF, with scale "RBF_scale"
% on some specified "centers".
% The RBF is the "RBF_type"-th in the list of functions described in
% "RBF.m" (see description below).
%
% The function is vectorial, i.e. the number of centers may be bigger than
% 1.
%
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
%
% RBF_scale: a positive scalar: all the RBFs have the same scale.
%            a vector of the same number of rows of centers: scale of the
%            RBF for each center.
%--------------------------------------------------------------------------
% Output:
%--------------------------------------------------------------------------
% moms: vector containing the moments, i.e. the "k"-th component contains
%      the moment of the RBF with center "centers(k,:)" and scale
%      "RBF_scale(k)", that is
%              integral(phi(norm(X-center(k,:))/RBF_scale(k)),P)
%       being "P" the polygonal domain, "phi" the RBF
%--------------------------------------------------------------------------
% Copyrights.
%--------------------------------------------------------------------------
%% Copyright (C) 2019-2020 Alvise Sommariva, Marco Vianello.
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
%% Date: November 2019-October 16, 2020.
%--------------------------------------------------------------------------

% ............................ key idea ...................................
% After the command "boundary(P)", we save in "X" and "Y" the absissae
% describing the polygonal domain.
%
% Each polygonal component describing the boundary is separated from the
% next one by a NaN.
%
% Example: if [xv1 yv1] tracks the outer boundary and [xv2 yv2] tracks the
% inner boundary (a hole), then
%           "X=[xv1; NaN; xv2]" and "Y=[yv1; NaN; yv2];"
%
% "Boundary" command gives clockwise orientation to the polygon component
% if it is not  a hole, counterclockwise otherwise; thus the contribution
% of the moments, thought respectively counterclockwise/clockwise must be
% changed of sign.
% .........................................................................

% .......................... Troubleshooting ................................

if nargin < 1, error('Domain not defined.'); end
if nargin < 2, error('Centers not defined.'); end
if nargin < 3, warning('RBF not declared. Using TPS.'); RBF_type=5; end
if nargin < 4,
    warning('RBF_scale not declared. Using RBF_scale=1.');
    RBF_scale=1;
end


% ................. polygon described as polyshape object .................

% note: the code below needs the polyshape subroutines (Matlab versions
% at least Matlab 9.3.0.713579).

[X,Y]=boundary(P);
inans=find(isnan(X) == 1);


% .............. analysing first components and their moments ..............

% Here we study each component of the domain.

i0=1;
for ii=1:length(inans)
    i1=inans(ii)-1;
    xvL=X(i0:i1); yvL=Y(i0:i1);
    momsL(:,ii)=RBF_moms(xvL,yvL,centers,RBF_type,RBF_scale);
    i0=i1+2;
end


% ................ analysing last component and its moments ...............

i1=length(X);
xvL=(X(i0:i1));
yvL=(Y(i0:i1));
iend=length(inans)+1;
momsL(:,iend)=RBF_moms(xvL,yvL,centers,RBF_type,RBF_scale);


% ........................ summing contributions ..........................

% change of sign due to polyshape description of the boundary
moms=-sum(momsL,2);










function moms=RBF_moms(xv,yv,centers,RBF_type,RBF_scale)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% Once it has been defined the polygon with vertices (xv,yv), this routine
% computes the moments of a RBF, with scale "RBF_scale".
% The RBF is the "RBF_type"-th in the list of functions described in
% "RBF.m".
%
% Note: The function is vectorial, i.e. the number of centers may be
% bigger than 1.
%--------------------------------------------------------------------------
% Input:
%--------------------------------------------------------------------------
% xv: abscissae of the vertices (column vector).
%
% yv: ordinates of the vertices (column vector).
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
% 9: phi=@(r) (max(0,(1-r))).^2;             % Wendland W0
% 10: phi=@(r) (35*r.^2+18*r+3).*(max(0,(1-r))).^6;         % Wendland W4
% 11: phi=@(r) (32*r.^3+25*r.^2+8*r+1).*(max(0,(1-r))).^8;  % Wendland W6
% 12: phi=@(r) (sqrt(2)/(3*sqrt(pi))*(3*r^2*log(r/(1+sqrt(1-r.^2)))+...
%            (2*r^2+1).*sqrt(1-r^2)))*max(0,(1-r));  % Missing Wendland
% 13 phi=@(r) exp(-r);                     % Matern beta_1=3/2.
% 14 phi=@(r) (1+r).*exp(-r);              % Matern beta_2=5/2.
%
% RBF_scale: RBF scale, possibly depending on the center.
%
%-------------------
% Important:
%-------------------
% Differently from classical codes on polygons, the first and last vertex
% may not be equal, i.e. it is not necessary that
%             [xv(1) yv(1)]=[xv(end),yv(end)].
%--------------------------------------------------------------------------
% Output:
%--------------------------------------------------------------------------
% moms: vector containing the moments, whose "k"-th component contains the
%     moment of the RBF with center "centers(k,:)" and scale "RBF_scale(k)".
%--------------------------------------------------------------------------
% Copyrights.
%--------------------------------------------------------------------------
%% Copyright (C) 2019-2020 Alvise Sommariva, Marco Vianello.
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
%% Date: November 13, 2019- 10 October 2020
%--------------------------------------------------------------------------

% Set as equal first and last vertex.
if norm([xv(1) yv(1)]-[xv(end) yv(end)]) > 0
    xv=[xv; xv(1)]; yv=[yv; yv(1)];
end

Nvertices=length(xv);
Ncenters=size(centers,1);

% .................... call primitive functions ...........................
[psi,primitive]=primitives_list(RBF_type);


% .................... compute side contributions .........................

% ............................ Key ideas ..................................
% The routine works as follows:
%
% 1. it fixes a side of the polygon with vertices "V1" and "V2";
% 2. it computes for each center the contribution to the moment given by
%    the side via evaluating the "primitive" integral on the segment;
% 3. it sums the contributions.
% .........................................................................

moms=zeros(Ncenters,1);

for k=1:Nvertices-1
    
    % ..... computing integrals over the k-th side .....
    
    % ... k-th side vertices ...
    V1=[xv(k) yv(k)];
    V2=[xv(k+1) yv(k+1)];
    
    
    % Knowing sides C_V1, C_V2, V1_v2 we determine the angles of C_V1_V2,
    % and next a parameter "c" fundamental to compute the desired integral.
    
    V1_C=centers-repmat(V1,Ncenters,1);
    dist_V1_C=sqrt(V1_C(:,1).^2+V1_C(:,2).^2);
    V2_C=centers-repmat(V2,Ncenters,1);
    dist_V2_C=sqrt(V2_C(:,1).^2+V2_C(:,2).^2);
    dist_V1_V2=norm(V1-V2);
    
    angle_V1_C_V2=acos((dist_V1_C.^2+dist_V2_C.^2-dist_V1_V2.^2)./...
        (2*dist_V1_C.*dist_V2_C));
    
    angle_C_V1_V2=acos((dist_V1_C.^2+dist_V1_V2.^2-dist_V2_C.^2)./...
        (2*dist_V1_C.*dist_V1_V2));
    
    c=dist_V1_C.*sin(angle_C_V1_V2);
    
    % ... integration extremal points ...
    t0=angle_C_V1_V2;
    t1=angle_V1_C_V2+angle_C_V1_V2;
    
    
    % Note: patch, needed in non convex polygons, due to possible change of
    % directions.  It may happen that the segment (V1-V2) is clockwise and
    % changes moments sign contribution by (V1-V2).
    
    for ii=1:Ncenters
        SL=cross([V1_C(ii,:) 0],[V2_C(ii,:) 0]);
        S(ii,1)=sign(SL(3));
    end
    
    switch RBF_type
        
        case {4,9,10,11}
            
            % ............... methods for compacted support RBF ...........
            
            c=c./RBF_scale;
            for ii=1:Ncenters
                cL=c(ii); T0=t0(ii); T1=t1(ii);
                value0=primitive_compact_support(T0,cL,psi,primitive);
                value1=primitive_compact_support(T1,cL,psi,primitive);
                momsL(ii,1)=value1-value0;
            end
            
        case 2
            
            % *** key idea ***
            % With regard to the gaussian RBF, a primitive is not available
            % and we try to compute the integrals
            %           int(psi(c/sin(t)),t0,t1)
            % via suitable Gauss-Legendre rules, whose order depend on "c",
            % "t0" and "t1".
            % This approach is faster then numerical adaptive integration.
            
            c=c./RBF_scale;
            momsL=RBFmoms_gauss(c,t0,t1,psi);
            
        case {12,13,14}
            
            % .............. methods for RBF without primitives ...........
            
            % ... using numerical integration ...
            c=c./RBF_scale;
            for ii=1:Ncenters
                cL=c(ii); T0=t0(ii); T1=t1(ii);
                psiL=@(t) psi(cL./sin(t))-psi(0);
                tol=10^(-20);
                momsL(ii,1)=integral(psiL,T0,T1,...
                    'AbsTol',tol,'RelTol',tol);
            end
            
        case {1,3}
            
            method_123=1;
            switch method_123
                case 1
                    % Separate "singular"/"regular" values of t0,t1,cL.
                    c=c./RBF_scale; mtol=10^(-3); Mtol=1; 
                    momsL=zeros(Ncenters,1);
                    iok=find(c > mtol & t0 > mtol & t1 > mtol ...
                        & c < Mtol & t0 < Mtol & t1 < Mtol);
                    iko=setdiff(1:Ncenters,iok);
                    t0_ok=t0(iok); t1_ok=t1(iok); c_ok=c(iok);
                    
                    % via primitive
                    momsL(iok)=(primitive(t1_ok,c_ok)-primitive(t0_ok,c_ok));
                    
                    % via adaptive integration
                    for ii=1:length(iko)
                        ikoL=iko(ii);
                        cL=c(ikoL); T0=t0(ikoL); T1=t1(ikoL);
                        psiL=@(t) psi(cL./sin(t))-psi(0);
                        tol=10^(-25);
                        momsL(ikoL,1)=integral(psiL,T0,T1,'AbsTol',tol,...
                            'RelTol',tol);
                        % momsLp(ikoL,1)=primitive(T1,cL)-primitive(T0,cL);
                        % AE_ko(ii)=abs((momsLp(ikoL,1)-momsL(ikoL,1)));
                        % RE_ko(ii)=abs((momsLp(ikoL,1)-momsL(ikoL,1))/momsL(ikoL,1));
                    end
                case 2
                    cL=c./RBF_scale;
                    momsL=(primitive(t1,cL)-primitive(t0,cL));
            end
            
        otherwise % case {5,6,7,8}
            % ............ methods for not compacted support RBF ..........
            
            cL=c./RBF_scale;
            momsL=(primitive(t1,cL)-primitive(t0,cL));
            
    end
    
    moms=moms+S.*momsL;
end

moms=(RBF_scale.^2).*moms;
moms=real(moms);









function value=primitive_compact_support(t,c,psi,primitive)


%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This subroutine computes the moment of a certain RBF with compact support
% whose primitive would be "primitive" if the support would not be compact.
% In particular "primitive" has derivative w.r.t. "t" equal to
% "psi(c/sin(t)".
%--------------------------------------------------------------------------
% INPUT
%--------------------------------------------------------------------------
% t: point in which evaluate the primitive of compact support RBF;
% c: parameter defined by center and polygon sides;
% psi: "inner" primitive of the RBF, w.r.t. the radius;
% primitive: primitive of the RBF if it would not be compactly supported.
%--------------------------------------------------------------------------

if c > 1
    value=psi(1)*(t-pi/2);
else
    asinc=asin(c);
    if (t >= asinc ) & (t <= pi-asinc )
        value=primitive(t,c)-primitive(pi/2,c);
    else
        if t < asinc
            value=primitive(asinc,c)-primitive(pi/2,c);
            value=value+psi(1)*(t-asinc);
        else
            %fprintf('\n \t norm(t-pi/2) >  asin(c), t >= pi/2');
            value=primitive(pi-asinc,c)-primitive(pi/2,c);
            value=value+psi(1)*(t-pi+asinc);
        end
    end
end










function [psi,primitive]=primitives_list(RBF_type)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% Computing RBF primitives. In the case of cubature with Gauss-Green
% theorem, if "phi" is the RBF, we compute:
%                       psi=int(r*phi(r),r)
%                 primitive=int(psi(c/sin(t)),t)
% for a certain parameter "c" depending on the center and on each side.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% % RBF_type: the choosen RBF is the "RBF_type"-th in the list of functions
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
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% psi: if "phi" is the RBF, then psi=int(r*phi(r),r)
% primitive: for a certain parameter "c" depending on the center and on
%         each side "primitive=int(psi(c/sin(t)),t)".
%       In the case of compact support RBF, it is the primitive if compact
%       support would not be taken into account.
%--------------------------------------------------------------------------
% Copyright.
%--------------------------------------------------------------------------
%% Copyright (C) 2019-2019 Alvise Sommariva, Marco Vianello.
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
%% Date: November 13, 2019.
%--------------------------------------------------------------------------

pert=@(t) eps*(1-abs(sign(t)));
primitive=@(t,c) zeros(size(t,1),1)+0*c;

switch RBF_type
    
    
    case 1
        
        % ....................... Multiquadric ............................
        
        phi=@(r) (1+r.*r).^(1/2);
        psi=@(r) (r.^2 + 1).^(3/2)/3;
        primitive=@(t,c) real(i.*log(sqrt(2.*c.^2-cos(2.*t)+1)+...
            i.*sqrt(2).*cos(t))-1./2.*sqrt(c.^2).*(c.^2+3).*...
            atanh((sqrt(2).*sqrt(c.^2).*cos(t))./...
            sqrt(2.*c.^2-cos(2.*t)+1))-1./2.*c.^2.*cot(t).*csc(t).*...
            sqrt(c.^2-1./2.*cos(2.*t)+1./2))/3-t/3;
        
        
    case 2
        
        % ......................... Gaussian ..............................
        
        %% available via numerics
        phi=@(r) exp(-r.*r);
        psi=@(r) -exp(-r.^2)/2;
        
        
    case 3
        
        % .................... Inverse Multiquadric .......................
        
        phi=@(r) (1+r.*r).^(-1/2);
        psi=@(r) (r.^2 + 1).^(1/2);
        
        term1=@(t,c) sin(t+pert(t)).*sqrt(2*c.^2.*(csc(t+pert(t))).^2+2);
        term2=@(t,c) sqrt(-c.^2).*atanh((sqrt(2)*sqrt(-c.^2).*cos(t))./...
            sqrt(-2*c.^2+cos(2*t)-1));
        term3=@(t,c) -log( sqrt(-2*c.^2+cos(2*t)-1)+ sqrt(2)*cos(t) );
        term4=@(t,c) sqrt(-2*c.^2+cos(2*t)-1);
        term5=@(t,c) psi(0)*t;
        primitive=@(t,c) -(term1(t,c).*(term2(t,c)+term3(t,c))./...
            term4(t,c))-term5(t,c);
        
        
    case 4
        
        % ..................... Wendland 2 ................................
        
        phi=@(r) (1+4.*r).*(max(0,(1-r)).^4);  %
        psi_pre=@(r) (r.^2.*(8.*r.^5 - 35.*r.^4 + 56.*r.^3 - ...
            35.*r.^2 + 7))/14;
        psi=@(r) (r <= 1).*(psi_pre(r))+(r > 1).*psi_pre(1);
        
        primitive=@(t,c) (7.*c.^6.*tan(t./2)-(2.*c.^7)./3+ ...
            (35.*c.^4.*tan(t./2).^3.*(5.*c.^2+4))./3- ...
            2.*c.^5.*tan(t./2).^2.*(3.*c.^2+14)- ...
            2.*c.^5.*tan(t./2).^4.*(15.*c.^2+112)+ ...
            14.*c.^2.*tan(t./2).^5.*(25.*c.^4+30.*c.^2-8))./...
            (448.*tan(t./2).^6)-(c.^6.*tan(t./2).^5)./64+ ...
            (c.^7.*tan(t./2).^6)./672-(c.^2.*tan(t./2).* ...
            (25.*c.^4+30.*c.^2-8))./32-(5.*c.^4.*tan(t./2).^3.* ...
            (5.*c.^2+4))./192+(c.^5.*tan(t./2).^4.*(3.*c.^2+14))./224+...
            (c.^5.*tan(t./2).^2.*(15.*c.^2+112))./224+ ...
            (c.^5.*log(tan(t./2)).*(5.*c.^2+42))./28;
        
        
    case 5
        % .................. Thin Plate Spline ............................
        
        pert=@(r) eps.*(1-abs(sign(r)));
        phi=@(r) r.*r.*log(r+pert(r));
        psi=@(r) (r.^4.*(log(r) - 1/4))/4;
        
        primitive=@(t,c) (c.^4.*t)/6+13/72*c.^4.*cot(t)+...
            7/144.*c.^4.*(cos(t)./(sin(t)+pert(t)).^3)+...
            1/12*c.^4.*(cos(2*t)-2).*(cos(t)./(sin(t)+pert(t)).^3).*...
            log(c.*(1./(sin(t)+pert(t))));
        
        
    case 6
        
        % ................... polyharmonic spline r^3 .....................
        
        phi=@(r) r.^3;
        psi=@(r) r.^5/5;
        
        primitive=@(t,c) (c.^5.*cos(t).*(3*cos(t).^2 - 5))./...
            (40*(cos(t).^2 - 1+pert(t)).^2) - (3*c.^5.*atanh(cos(t)))./40;
        
    case 7
        
        % ................... polyharmonic spline r^5 .....................
        
        phi=@(r) r.^5;
        psi=@(r) r.^7/7;
        
        primitive=@(t,c)(c.^7.*cos(t).*(15*cos(t).^4-40*cos(t).^2+33))./...
            (336*(cos(t).^2 - 1+pert(t)).^3) - (5*c.^7.*atanh(cos(t)))/112;
        
    case 8
        
        % ................... polyharmonic spline r^7 .....................
        
        phi=@(r) r.^7;
        psi=@(r) r.^9/9;
        
        primitive=@(t,c) (c.^9.*cos(t).*(511.*cos(t).^2-385.*cos(t).^4+...
            105.*cos(t).^6 - 279))./(3456.*(cos(t).^2 - 1).^4) - ...
            (35.*c.^9.*atanh(cos(t)))./1152;
        
    case 9
        % ....................... Wendland W0 .............................
        phi=@(r) max(0,(1-r).^2);
        psi_pre=@(r) (r.^2.*(3.*r.^2 - 8.*r + 6))/12;
        psi=@(r) (r <= 1).*(psi_pre(r))+(r > 1).*psi_pre(1);
        
        primitive=@(t,c) tan(t./2).*((3*c.^4)/32 + c.^2/4) - ...
            (tan(t./2).^2.*(9*c.^4 + 24.*c.^2) - ...
            8.*c.^3.*tan(t./2) + c.^4)./(96.*tan(t./2).^3) - ...
            (c.^3.*tan(t./2).^2)/12 + (c.^4.*tan(t./2).^3)/96 - ...
            (c.^3.*log(tan(t./2)))./3;
        
    case 10
        
        % ....................... Wendland W4 .............................
        
        phi=@(r) (35.*r.^2+18.*r+3).*max(0,(1-r).^6);
        psi_pre=@(r) ((r.^2.*(21.*r.^8 - 128.*r.^7 + 315.*r.^6 - ...
            384.*r.^5 + 210.*r.^4 - 42.*r.^2 + 9))/6);
        psi=@(r) (r <= 1).*(psi_pre(r))+(r > 1).*(r-1).*psi_pre(1);
        
        primitive=@(t,c) -64./45.*c.^10.*cot(t)-7./18.*c.^10.*cot(t).* ...
            (csc(t)).^8-4./9.*c.^10.*cot(t).*(csc(t)).^6- ...
            8./15.*c.^10.*cot(t).*(csc(t)).^4-32./45.*c.^10.*cot(t).* ...
            (csc(t)).^2+1./96.*c.^9.*(csc(t./2)).^8+5./72.*c.^9.* ...
            (csc(t./2)).^6+5./16.*c.^9.*(csc(t./2)).^4+35./24.* ...
            c.^9.*(csc(t./2)).^2-1./96.*c.^9.*(sec(t./2)).^8- ...
            5./72.*c.^9.*(sec(t./2)).^6-5./16.*c.^9.*(sec(t./2)).^4- ...
            35./24.*c.^9.*(sec(t./2)).^2-35./6.*c.^9.*log(sin(t./2))+ ...
            35./6.*c.^9.*log(cos(t./2))-24.*c.^8.*cot(t)- ...
            15./2.*c.^8.*cot(t).*(csc(t)).^6-9.*c.^8.*cot(t).* ...
            (csc(t)).^4-12.*c.^8.*cot(t).*(csc(t)).^2+1./6.* ...
            c.^7.*(csc(t./2)).^6+c.^7.*(csc(t./2)).^4+5.*c.^7.* ...
            (csc(t./2)).^2-1./6.*c.^7.*(sec(t./2)).^6-c.^7.* ...
            (sec(t./2)).^4-5.*c.^7.*(sec(t./2)).^2-20.*c.^7.* ...
            log(sin(t./2))+20.*c.^7.*log(cos(t./2))-56./3.*c.^6.* ...
            cot(t)-7.*c.^6.*cot(t).*(csc(t)).^4-28./3.*c.^6.*cot(t).* ...
            (csc(t)).^2+14./3.*c.^4.*cot(t)+7./3.*c.^4.*cot(t).* ...
            (csc(t)).^2-3./2.*c.^2.*cot(t);
        
    case 11
        
        % ....................... Wendland W6 .............................
        
        phi=@(r) (32.*r.^3+25.*r.^2+8.*r+1).*max(0,(1-r).^8);
        psi_pre=@(r) (32.*r.^13)/13 - (77.*r.^12)/4 + 64.*r.^11 - ...
            (231.*r.^10)/2 + (352.*r.^9)/3 - (231.*r.^8)/4 + ...
            11.*r.^6 - (11.*r.^4)/4 + r.^2/2;
        psi=@(r) (r <= 1).*(psi_pre(r))+(r > 1).*(r-1).*psi_pre(1);
        
        primitive=@(t,c) -(c.^13.*(csc(t./2)).^12)./19968-...
            (7.*c.^13.*(csc(t./2)).^10)./16640-...
            (7.*c.^13.*(csc(t./2)).^8)./3328-7./832.*c.^13.*...
            (csc(t./2)).^6-(105.*c.^13.*(csc(t./2)).^4)./...
            3328-(231.*c.^13.*(csc(t./2)).^2)./1664+(c.^13.*...
            (sec(t./2)).^12)./19968+(7.*c.^13.*(sec(t./2)).^10)./...
            16640+(7.*c.^13.*(sec(t./2)).^8)./3328+7./832.*c.^13.* ...
            (sec(t./2)).^6+(105.*c.^13.*(sec(t./2)).^4)./3328+ ...
            (231.*c.^13.*(sec(t./2)).^2)./1664+231./416.*c.^13.* ...
            log(sin(t./2))-231./416.*c.^13.*log(cos(t./2))+64./9.* ...
            c.^12.*cot(t)+7./4.*c.^12.*cot(t).*(csc(t)).^10+35./18.* ...
            c.^12.*cot(t).*(csc(t)).^8+20./9.*c.^12.*cot(t).* ...
            (csc(t)).^6+8./3.*c.^12.*cot(t).*(csc(t)).^4+ ...
            32./9.*c.^12.*cot(t).*(csc(t)).^2-1./160.*c.^11.* ...
            (csc(t./2)).^10-3./64.*c.^11.*(csc(t./2)).^8-7./32.* ...
            c.^11.*(csc(t./2)).^6-7./8.*c.^11.*(csc(t./2)).^4-63./16.* ...
            c.^11.*(csc(t./2)).^2+1./160.*c.^11.*(sec(t./2)).^10+ ...
            3./64.*c.^11.*(sec(t./2)).^8+7./32.*c.^11.*(sec(t./2)).^6+ ...
            7./8.*c.^11.*(sec(t./2)).^4+63./16.*c.^11.*(sec(t./2)).^2+ ...
            63./4.*c.^11.*log(sin(t./2))-63./4.*c.^11.*log(cos(t./2))+ ...
            704./15.*c.^10.*cot(t)+77./6.*c.^10.*cot(t).*(csc(t)).^8+ ...
            44./3.*c.^10.*cot(t).*(csc(t)).^6+88./5.*c.^10.*cot(t).* ...
            (csc(t)).^4+352./15.*c.^10.*cot(t).*(csc(t)).^2- ...
            11./192.*c.^9.*(csc(t./2)).^8-55./144.*c.^9.*...
            (csc(t./2)).^6-55./32.*c.^9.*(csc(t./2)).^4- ...
            385./48.*c.^9.*(csc(t./2)).^2+11./192.*c.^9.* ...
            (sec(t./2)).^8+55./144.*c.^9.*(sec(t./2)).^6+55./32.* ...
            c.^9.*(sec(t./2)).^4+385./48.*c.^9.*(sec(t./2)).^2+ ...
            385./12.*c.^9.*log(sin(t./2))-385./12.*c.^9.*log(cos(t./2))+...
            132./5.*c.^8.*cot(t)+33./4.*c.^8.*cot(t).*(csc(t)).^6+99./ ...
            10.*c.^8.*cot(t).*(csc(t)).^4+66./5.*c.^8.*cot(t).* ...
            (csc(t)).^2-88./15.*c.^6.*cot(t)-11./5.*c.^6.*cot(t).* ...
            (csc(t)).^4-44./15.*c.^6.*cot(t).*(csc(t)).^2+11./6.* ...
            c.^4.*cot(t)+11./12.*c.^4.*cot(t).*(csc(t)).^2-1/2.* ...
            c.^2.*cot(t);
        
    case 12 %% Not working.
        
        % ..................... Missing Wendland ..........................
        
        pert=@(r) eps*(1-abs(sign(r)));
        
        % 12: phi=@(r) (sqrt(2)/(3*sqrt(pi))*(3*r^2*log(r/(1+sqrt(1-r.^2)))+...
        %            (2*r^2+1).*sqrt(1-r^2)))*max(0,(1-r));
        
        phi=@(r) (...
            sqrt(2)/(3.*sqrt(pi)).*...
            (...
            3.*r.^2.*log( r./(1+sqrt(1-r.^2)) + pert(r) )+...
            (2.*r.^2+1).*sqrt(1-r.^2))...
            ).*max(0,(1-r));
        
        psi_pre=@(r) 1/3*sqrt(2/pi)*(1/20*sqrt(1-r.^2).*...
            (8*r.^4+9*r.^2-2)+3/4*r.^4.*log(pert(r)+r./(sqrt(1-r.^2)+1)));
        
        psi_pre=@(r) psi_pre(r);
        
        %         psi_pre=@(r) (1 - r.^2).^(1/2).*((2.*r.^4)/5 + ...
        %             (9.*r.^2)/20 - 1/10) + (3.*r.^4.* ...
        %             log((r+pert(r))/((1 - r.^2).^(1/2) + 1)))/4;
        %
        %         psi_pre=@(r) sqrt(2)/(3.*sqrt(pi))* psi_pre(r);
        
        psi=@(r) (r <= 1).*(psi_pre(r))+(r > 1).*psi_pre(1);
        
        primitive_pre=@(t,c) (1/40)*(-5*c.^2.*cot(t).*sqrt((2.*c.^2+ ...
            cos(2.*t)-1)./(cos(2.*t)-1))+1./2.*c.^2.*cot(t).*...
            (csc(t)).^2.*((6.*c.^2+7).*cos(2.*t)-7.*(2.*c.^2+1)).*...
            sqrt(1-c.^2.*(csc(t)).^2)+(5.*sqrt(c.^2).*csc(t).* ...
            sqrt(2.*c.^2+cos(2.*t)-1).*((3.*c.^2+1).* ...
            atanh((sqrt(2).*sqrt(c.^2).*cos(t))./...
            sqrt(2.*c.^2+cos(2.*t)-1))-4.*(c.^2).^(3./2).* ...
            atanh((sqrt(2).*cos(t))./sqrt(2.*c.^2+cos(2.*t)-1))))./...
            sqrt((csc(t)).^2.*(-(2.*c.^2+cos(2.*t)-1)))+...
            10.*c.^4.*(cos(2.*t)-2).*cot(t).*(csc(t)).^2.*...
            log((c.*csc(t))./(sqrt((2.*c.^2+cos(2.*t)-1)./...
            (cos(2.*t)-1))+1))-(sin(t).*sqrt(2-2.*c.^2.*(csc(t)).^2).*...
            (4.*sqrt(c.^2).*log(sqrt(2.*c.^2+cos(2.*t)-1)+...
            sqrt(2).*cos(t))+(6.*c.^4+5.*c.^2-15).*c.^2.*...
            atanh((sqrt(2).*sqrt(c.^2).*cos(t))./ ...
            sqrt(2.*c.^2+cos(2.*t)-1))))./(sqrt(c.^2).*sqrt(2.*c.^2+...
            cos(2.*t)-1)));
        
        pertpi=@(t) -eps*(1-abs(sign(t-pi/2)));
        primitive=@(t,c) primitive_pre(t+pertpi(t),c)-...
            (-2/(30*sqrt(2*pi)))*t;
        
    case 13
        
        % ................ Matern, beta_1=(d+1)/2, d=2; ...................
        
        phi=@(r) exp(-r);
        psi=@(r) -exp(-r).*(r + 1)+1;
        
    case 14
        
        % ................ Matern beta_2=(d+3)/2, where d=2 ...............
        
        phi=@(r) (1+r).*exp(-r);
        psi=@(r) -exp(-r).*(r.^2 + 3*r + 3)-1;
        
        
    otherwise
        error('RBF type not implemented')
end










function phi =  RBF (RBF_type)

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
%--------------------------------------------------------------------------
% Output:
%--------------------------------------------------------------------------
% phi: RBF function.
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
        phi=@(r) (1+r.*r).^(1/2);            % Multiquadric
    case 2
        phi=@(r) exp(-r.*r);                  % Gaussian
    case 3
        phi=@(r) (1+r.*r).^(-1/2);            % Inverse Multiquadric
    case 4
        phi=@(r) (1+4*r).*(max(0,(1-r))).^4;  % Wendland 2
    case 5
        pert=@(r) eps*(1-abs(sign(r)));
        phi=@(r) r.*r.*log(r+pert(r)); % TPS
    case 6
        phi=@(r) r.^3;                        % polyharmonic spline
    case 7
        phi=@(r) r.^5;                        % polyharmonic spline
    case 8
        phi=@(r) r.^7;                        % polyharmonic spline
    case 9
        phi=@(r) (max(0,(1-r))).^2;             % Wendland W0
    case 10
        phi=@(r) (35*r.^2+18*r+3).*(max(0,(1-r))).^6;  % Wendland W4
    case 11
        phi=@(r) (32*r.^3+25*r.^2+8*r+1).*(max(0,(1-r))).^8;  % Wendland W6
    case 12                                    % Missing Wendland
        pert=@(r) eps*(1-abs(sign(r)));
        phi=@(r) (sqrt(2)/(3*sqrt(pi))*...
            (3*r.^2.*log( r./(1+sqrt(1-r.^2)) + pert(r) )+...
            (2*r.^2+1).*sqrt(1-r.^2))).*max(0,(1-r));
    case 13                             % Matern beta_1=(d+1)/2, where d=2.
        phi=@(r) exp(-r);
    case 14                             % Matern beta_2=(d+3)/2, where d=2.
        phi=@(r) (1+r).*exp(-r);
        
    otherwise
        error('RBF type not implemented')
end











function moms=RBFmoms_gauss(c,t0,t1,psi)

%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% In this procedure we compute the integral
%
%               integral(psi(c/sin(t)),t0,t1)
%
% via Gauss-Legendre rules, whose degree depends on the integrand.
%
% We have fixed some Gauss-Legendre rules, with degree of precision
%
%                    N=20,100,250,500,1400
%
% and use them coherently with the integrand.
%
% The degree N has been choosen by determining the length M of a chebfun
% function approximating the integrand above at machine precision, and then
% setting N >= M.
%
% If the degree was over 1400, we computed the integral by means of Matlab
% built-in adaptive routine.
%
% The so obtained values of the integrand, allows to compute quickly the
% RBF moments of Gaussian RBF over a polygon.
%
% In particular, with regard to the gaussian RBF, a primitive is not
% available and we try to compute the integrals
%           int(psi(c/sin(t)),t0,t1)
% via suitable Gauss-Legendre rules, whose order depend on "c",
% "t0" and "t1".  This approach is faster then numerical adaptive
% integration.
%--------------------------------------------------------------------------
% INPUT:
%
% c,t0,t1,psi: parameters needed to determine the integrals to be computed,
%   i.e. integral(psi(c/sin(t)),t0,t1).
%--------------------------------------------------------------------------
% OUTPUT:
%
% moms: approximation of "integral(psi(c/sin(t)),t0,t1)" to achieve moments
%  of gaussians RBF.
%--------------------------------------------------------------------------
% Copyrights.
%--------------------------------------------------------------------------
%% Copyright (C) 2019-2020 Alvise Sommariva, Marco Vianello.
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
%% Date: October 2019 - October 10, 2020.
%--------------------------------------------------------------------------

[xXS,wXS]=gauss_legendre_10;
[xS,wS]=gauss_legendre_50;
[xSM,wSM]=gauss_legendre_125;
[xM,wM]=gauss_legendre_250;
[xL,wL]=gauss_legendre_700;

for ii=1:length(t0)
    cL=c(ii); T0=t0(ii); T1=t1(ii);
    
    psiL=@(t) psi(cL./sin(t))-psi(0);
    
    moms(ii,1)=compute_integral(cL,T0,T1,psiL,xXS,wXS,xS,wS,xSM,wSM,...
        xM,wM,xL,wL);
    
end










function val=compute_integral(cL,T0,T1,psiL,xXS,wXS,xS,wS,xSM,wSM,...
    xM,wM,xL,wL)


if cL <= 10^(-8)
    if T0 > 1/5 & T1 < pi-1/5
        % fprintf('\n \t XSMALL 1');
        val=gaussian_cubature(T0,T1,xXS,wXS,psiL);
    else
        if T0 > 1/1000 & T1 < pi-1/1000
            % fprintf('\n \t SMALL 1');
            val=gaussian_cubature(T0,T1,xS,wS,psiL);
        else
            % fprintf('\n \t NUMERICAL 1');
            tol=10^(-20); val=integral(psiL,T0,T1,'AbsTol',tol,'RelTol',tol);
        end
        return;
    end
end


if cL > 10^(-8) & cL <= 10^(-2)
    if T0 > 1/10 & T1 < pi-1/10
        % fprintf('\n \t SMALL 2');
        val=gaussian_cubature(T0,T1,xS,wS,psiL);
    else
        if T0 > 1/100 & T1 < pi-1/100
            % fprintf('\n \t MEDIUM 2');
            val=gaussian_cubature(T0,T1,xM,wM,psiL);
        else
            if T0 > 1/1000 & T1 < pi-1/1000
                % fprintf('\n \t LARGE 2');
                val=gaussian_cubature(T0,T1,xL,wL,psiL);
            else
                % fprintf('\n \t NUMERICAL 2');
                tol=10^(-20);
                val=integral(psiL,T0,T1,'AbsTol',tol,'RelTol',tol);
            end
        end
    end
    return;
end


if cL > 10^(-2) & cL <= 10^(-1)
    if T0 > 1/100 & T1 < pi-1/100
        % fprintf('\n \t MEDIUM 3 SM');
        val=gaussian_cubature(T0,T1,xM,wM,psiL);
    else
        if T0 > 1/1000 & T1 < pi-1/1000
            % fprintf('\n \t LARGE 3 SM');
            val=gaussian_cubature(T0,T1,xL,wL,psiL);
        else
            % fprintf('\n \t NUMERICAL 3 SM');
            tol=10^(-20);
            val=integral(psiL,T0,T1,'AbsTol',tol,'RelTol',tol);
        end
        
    end
    return;
end



if cL > 10^(-1) & cL <= 10^(-0.5)
    if T0 > 1/10 & T1 < pi-1/10
        % fprintf('\n \t MEDIUM 3 MID safe');
        val=gaussian_cubature(T0,T1,xSM,wSM,psiL);
    else
        if T0 > 1/1000 & T1 < pi-1/1000
            % fprintf('\n \t MEDIUM 3 MID unsafe');
            val=gaussian_cubature(T0,T1,xM,wM,psiL);
        else
            % fprintf('\n \t NUMERICAL 3 MID');
            tol=10^(-20);
            val=integral(psiL,T0,T1,'AbsTol',tol,'RelTol',tol);
        end
        
    end
    return;
end





if cL > 10^(-0.5) & cL <= 26
    if T0 > 1/10000 & T1 < pi-1/10000
        % fprintf('\n \t MEDIUM 4 safe');
        val=gaussian_cubature(T0,T1,xM,wM,psiL);
    else
        % fprintf('\n \t NUMERICAL 4');
        tol=10^(-20);
        val=integral(psiL,T0,T1,'AbsTol',tol,'RelTol',tol);
    end
    return;
end



if cL > 26 & cL <= 28
    % fprintf('\n \t NUMERICAL 6');
    tol=10^(-20);
    val=integral(psiL,T0,T1,'AbsTol',tol,'RelTol',tol);
    return;
end


if cL > 28
    % fprintf('\n \t XS 7 safe');
    val=gaussian_cubature(T0,T1,xXS,wXS,psiL);
    return;
end










function val=gaussian_cubature(T0,T1,xG,wG,psiL)

xLOC=(T0+T1)/2+(T1-T0)*xG/2;
wLOC=(T1-T0)*wG/2;
fxLOC=feval(psiL,xLOC);
val=wLOC'*fxLOC;










function [X,W]=gauss_legendre_10

xw=[ -9.7390652851717174343093575e-01 6.6671344308688207380697577e-02
    -8.6506336668898453634568568e-01 1.4945134915058036484403203e-01
    -6.7940956829902443558921732e-01 2.1908636251598198607659640e-01
    -4.3339539412924721339948064e-01 2.6926671930999640514059479e-01
    -1.4887433898163118795032744e-01 2.9552422471475281451347428e-01
    1.4887433898163118795032744e-01 2.9552422471475281451347428e-01
    4.3339539412924721339948064e-01 2.6926671930999640514059479e-01
    6.7940956829902443558921732e-01 2.1908636251598198607659640e-01
    8.6506336668898453634568568e-01 1.4945134915058036484403203e-01
    9.7390652851717174343093575e-01 6.6671344308688207380697577e-02];

X=xw(:,1); W=xw(:,2);










function [X,W]=gauss_legendre_50

xw=[ -9.9886640442007101903243438e-01 2.9086225531551250338135883e-03
    -9.9403196943209071179126113e-01 6.7597991957454696757001678e-03
    -9.8535408404800584047933398e-01 1.0590548383650987690485223e-02
    -9.7286438510669204227099272e-01 1.4380822761485645075452133e-02
    -9.5661095524280792545823715e-01 1.8115560713489461258651758e-02
    -9.3665661894487795002817165e-01 2.1780243170124811286081368e-02
    -9.1307855665579185089342218e-01 2.5360673570012402106010896e-02
    -8.8596797952361305839019678e-01 2.8842993580535280367938000e-02
    -8.5542976942994608524628575e-01 3.2213728223578062814791423e-02
    -8.2158207085933598889937457e-01 3.5459835615146158283028655e-02
    -7.8455583290039931920745175e-01 3.8568756612587656862345398e-02
    -7.4449430222606849394395567e-01 4.1528463090147654801498334e-02
    -7.0155246870682230753146769e-01 4.4327504338803294658966081e-02
    -6.5589646568543935600814621e-01 4.6955051303948461272064208e-02
    -6.0770292718495022565861063e-01 4.9400938449466295920853298e-02
    -5.5715830451465009343081647e-01 5.1655703069581088149320180e-02
    -5.0445814490746421210332073e-01 5.3710621888996189221554545e-02
    -4.4980633497403876841502779e-01 5.5557744806212526478272906e-02
    -3.9341431189756514985589320e-01 5.7189925647728345747822232e-02
    -3.3550024541943734845972358e-01 5.8600849813222388728917167e-02
    -2.7628819377953200975284176e-01 5.9785058704265474360806110e-02
    -2.1600723687604175826670883e-01 6.0737970841770239083245997e-02
    -1.5489058999814589445698232e-01 6.1455899590316706571080374e-02
    -9.3174701560086142793082331e-02 6.1936067420683256490310242e-02
    -3.1098338327188876362150438e-02 6.2176616655347280437915458e-02
    3.1098338327188876362150438e-02 6.2176616655347280437915458e-02
    9.3174701560086142793082331e-02 6.1936067420683256490310242e-02
    1.5489058999814589445698232e-01 6.1455899590316706571080374e-02
    2.1600723687604175826670883e-01 6.0737970841770239083245997e-02
    2.7628819377953200975284176e-01 5.9785058704265474360806110e-02
    3.3550024541943734845972358e-01 5.8600849813222388728917167e-02
    3.9341431189756514985589320e-01 5.7189925647728345747822232e-02
    4.4980633497403876841502779e-01 5.5557744806212526478272906e-02
    5.0445814490746421210332073e-01 5.3710621888996189221554545e-02
    5.5715830451465009343081647e-01 5.1655703069581088149320180e-02
    6.0770292718495022565861063e-01 4.9400938449466295920853298e-02
    6.5589646568543935600814621e-01 4.6955051303948461272064208e-02
    7.0155246870682230753146769e-01 4.4327504338803294658966081e-02
    7.4449430222606849394395567e-01 4.1528463090147654801498334e-02
    7.8455583290039931920745175e-01 3.8568756612587656862345398e-02
    8.2158207085933598889937457e-01 3.5459835615146158283028655e-02
    8.5542976942994608524628575e-01 3.2213728223578062814791423e-02
    8.8596797952361305839019678e-01 2.8842993580535280367938000e-02
    9.1307855665579185089342218e-01 2.5360673570012402106010896e-02
    9.3665661894487795002817165e-01 2.1780243170124811286081368e-02
    9.5661095524280792545823715e-01 1.8115560713489461258651758e-02
    9.7286438510669204227099272e-01 1.4380822761485645075452133e-02
    9.8535408404800584047933398e-01 1.0590548383650987690485223e-02
    9.9403196943209071179126113e-01 6.7597991957454696757001678e-03
    9.9886640442007101903243438e-01 2.9086225531551250338135883e-03];

X=xw(:,1); W=xw(:,2);










function [X,W]=gauss_legendre_125

xw=[-9.9981641629776851765143419e-01 4.7112064849352385670355758e-04
    -9.9903283475297643967394379e-01 1.0963927055869938497617566e-03
    -9.9762362681096161676208567e-01 1.7219042078884714240538667e-03
    -9.9558936225216276838523299e-01 2.3464172536737898204506347e-03
    -9.9293127736660247162348014e-01 2.9694762015580532070468944e-03
    -9.8965102853875064337074718e-01 3.5906793228245121804564910e-03
    -9.8575066811248723830374274e-01 4.2096343701267004144828121e-03
    -9.8123263879985989088794440e-01 4.8259524859394437382165144e-03
    -9.7609977100610889610976528e-01 5.4392470750795895995111096e-03
    -9.7035528067547549557758657e-01 6.0491336609403379906413356e-03
    -9.6400276712800458955854310e-01 6.6552299959851358643336816e-03
    -9.5704621073949713849771115e-01 7.2571562508376178876612350e-03
    -9.4948997041646066019637829e-01 7.8545352306447387136234539e-03
    -9.4133878084891553505997308e-01 8.4469926013376224471773668e-03
    -9.3259774953497076577235703e-01 9.0341571191539012802840247e-03
    -9.2327235357550385685954097e-01 9.6156608605994613181433550e-03
    -9.1336843623924890422927092e-01 1.0191139451515405026094108e-02
    -9.0289220329957464716841287e-01 1.0760232294543732936564773e-02
    -8.9185021914478124216429933e-01 1.1322582794568035804982920e-02
    -8.8024940266409212874521018e-01 1.1877838581841214035672571e-02
    -8.6809702291176782384951593e-01 1.2425651732580315217413514e-02
    -8.5540069455196965364507378e-01 1.2965678986845308173769808e-02
    -8.4216837308716985255330201e-01 1.3497581963538951566050628e-02
    -8.2840834987306577463783697e-01 1.4021027372377320094343212e-02
    -8.1412924692309651675259374e-01 1.4535687222688631684008875e-02
    -7.9934001150580491490416080e-01 1.5041239028903245766866092e-02
    -7.8404991053841821546654955e-01 1.5537366012602575604528710e-02
    -7.6826852478015528191690464e-01 1.6023757300997306524115160e-02
    -7.5200574282889376398486547e-01 1.6500108121709166492108167e-02
    -7.3527175492495389086400337e-01 1.6966119993733048965101062e-02
    -7.1807704656588111635784344e-01 1.7421500914458813713547869e-02
    -7.0043239193622675031747349e-01 1.7865965542635340451704806e-02
    -6.8234884715644028574388358e-01 1.8299235377161517795974177e-02
    -6.6383774335510559172490730e-01 1.8721038931592100945655588e-02
    -6.4491067956885772538555557e-01 1.9131111904248723765142870e-02
    -6.2557951547443257922509474e-01 1.9529197343829424876604506e-02
    -6.0585636395740061210801741e-01 1.9915045810412725785232269e-02
    -5.8575358352224038416267149e-01 2.0288415531755431969740400e-02
    -5.6528377054850909022576388e-01 2.0649072554786002814397605e-02
    -5.4445975139796287667337538e-01 2.0996790892198663458501073e-02
    -5.2329457437756921045490799e-01 2.1331352664056231482891945e-02
    -5.0180150156344915934880646e-01 2.1652548234313024472230680e-02
    -4.7999400049087537212244570e-01 2.1960176342172229441151998e-02
    -4.5788573571552554364316734e-01 2.2254044228195393739788699e-02
    -4.3549056025128685121217131e-01 2.2533967755085031192674450e-02
    -4.1282250688997107479494275e-01 2.2799771523064654016321740e-02
    -3.8989577940838221481456571e-01 2.3051288979783902954867614e-02
    -3.6672474366824636682338223e-01 2.3288362524679927928472623e-02
    -3.4332391861457478565711199e-01 2.3510843607729623588875612e-02
    -3.1970796717811067466641362e-01 2.3718592822530807501246741e-02
    -2.9589168708755309022961910e-01 2.3911479993653991793500779e-02
    -2.7189000159731391281781043e-01 2.4089384258210177341963387e-02
    -2.4771795013662506468321567e-01 2.4252194141583249820115498e-02
    -2.2339067888584807075602612e-01 2.4399807627279869459702155e-02
    -1.9892343128589207168488429e-01 2.4532132220852823772938578e-02
    -1.7433153848669447061325855e-01 2.4649085007857923779184262e-02
    -1.4963040974073463229565562e-01 2.4750592705808066584793892e-02
    -1.2483552274761665346058948e-01 2.4836591710091889678713173e-02
    -9.9962413955758125383432855e-02 2.4907028133828425231488524e-02
    -7.5026668827269674122426579e-02 2.4961857841632465598857848e-02
    -5.0043912072120007306086364e-02 2.5001046477269750489824673e-02
    -2.5029797857712903635940549e-02 2.5024569485184552236622935e-02
    0.0000000000000000000000000e+00 2.5032412125886149140141512e-02
    2.5029797857712903635940549e-02 2.5024569485184552236622935e-02
    5.0043912072120007306086364e-02 2.5001046477269750489824673e-02
    7.5026668827269674122426579e-02 2.4961857841632465598857848e-02
    9.9962413955758125383432855e-02 2.4907028133828425231488524e-02
    1.2483552274761665346058948e-01 2.4836591710091889678713173e-02
    1.4963040974073463229565562e-01 2.4750592705808066584793892e-02
    1.7433153848669447061325855e-01 2.4649085007857923779184262e-02
    1.9892343128589207168488429e-01 2.4532132220852823772938578e-02
    2.2339067888584807075602612e-01 2.4399807627279869459702155e-02
    2.4771795013662506468321567e-01 2.4252194141583249820115498e-02
    2.7189000159731391281781043e-01 2.4089384258210177341963387e-02
    2.9589168708755309022961910e-01 2.3911479993653991793500779e-02
    3.1970796717811067466641362e-01 2.3718592822530807501246741e-02
    3.4332391861457478565711199e-01 2.3510843607729623588875612e-02
    3.6672474366824636682338223e-01 2.3288362524679927928472623e-02
    3.8989577940838221481456571e-01 2.3051288979783902954867614e-02
    4.1282250688997107479494275e-01 2.2799771523064654016321740e-02
    4.3549056025128685121217131e-01 2.2533967755085031192674450e-02
    4.5788573571552554364316734e-01 2.2254044228195393739788699e-02
    4.7999400049087537212244570e-01 2.1960176342172229441151998e-02
    5.0180150156344915934880646e-01 2.1652548234313024472230680e-02
    5.2329457437756921045490799e-01 2.1331352664056231482891945e-02
    5.4445975139796287667337538e-01 2.0996790892198663458501073e-02
    5.6528377054850909022576388e-01 2.0649072554786002814397605e-02
    5.8575358352224038416267149e-01 2.0288415531755431969740400e-02
    6.0585636395740061210801741e-01 1.9915045810412725785232269e-02
    6.2557951547443257922509474e-01 1.9529197343829424876604506e-02
    6.4491067956885772538555557e-01 1.9131111904248723765142870e-02
    6.6383774335510559172490730e-01 1.8721038931592100945655588e-02
    6.8234884715644028574388358e-01 1.8299235377161517795974177e-02
    7.0043239193622675031747349e-01 1.7865965542635340451704806e-02
    7.1807704656588111635784344e-01 1.7421500914458813713547869e-02
    7.3527175492495389086400337e-01 1.6966119993733048965101062e-02
    7.5200574282889376398486547e-01 1.6500108121709166492108167e-02
    7.6826852478015528191690464e-01 1.6023757300997306524115160e-02
    7.8404991053841821546654955e-01 1.5537366012602575604528710e-02
    7.9934001150580491490416080e-01 1.5041239028903245766866092e-02
    8.1412924692309651675259374e-01 1.4535687222688631684008875e-02
    8.2840834987306577463783697e-01 1.4021027372377320094343212e-02
    8.4216837308716985255330201e-01 1.3497581963538951566050628e-02
    8.5540069455196965364507378e-01 1.2965678986845308173769808e-02
    8.6809702291176782384951593e-01 1.2425651732580315217413514e-02
    8.8024940266409212874521018e-01 1.1877838581841214035672571e-02
    8.9185021914478124216429933e-01 1.1322582794568035804982920e-02
    9.0289220329957464716841287e-01 1.0760232294543732936564773e-02
    9.1336843623924890422927092e-01 1.0191139451515405026094108e-02
    9.2327235357550385685954097e-01 9.6156608605994613181433550e-03
    9.3259774953497076577235703e-01 9.0341571191539012802840247e-03
    9.4133878084891553505997308e-01 8.4469926013376224471773668e-03
    9.4948997041646066019637829e-01 7.8545352306447387136234539e-03
    9.5704621073949713849771115e-01 7.2571562508376178876612350e-03
    9.6400276712800458955854310e-01 6.6552299959851358643336816e-03
    9.7035528067547549557758657e-01 6.0491336609403379906413356e-03
    9.7609977100610889610976528e-01 5.4392470750795895995111096e-03
    9.8123263879985989088794440e-01 4.8259524859394437382165144e-03
    9.8575066811248723830374274e-01 4.2096343701267004144828121e-03
    9.8965102853875064337074718e-01 3.5906793228245121804564910e-03
    9.9293127736660247162348014e-01 2.9694762015580532070468944e-03
    9.9558936225216276838523299e-01 2.3464172536737898204506347e-03
    9.9762362681096161676208567e-01 1.7219042078884714240538667e-03
    9.9903283475297643967394379e-01 1.0963927055869938497617566e-03
    9.9981641629776851765143419e-01 4.7112064849352385670355758e-04];

X=xw(:,1); W=xw(:,2);










function [X,W]=gauss_legendre_250

xw=[-9.9995391943564204684236074e-01 1.1825670068171721075345887e-04
    -9.9975721221158908580406433e-01 2.7526103865831749294840192e-04
    -9.9940335329301610567398484e-01 4.3245460060038212139060798e-04
    -9.9889231972585434959910344e-01 5.8960034219324988433164059e-04
    -9.9822418225691300630586511e-01 7.4665740556624201132440710e-04
    -9.9739904367225573622590673e-01 9.0359824806269411524850543e-04
    -9.9641703298482942052771705e-01 1.0603974326420877692667144e-03
    -9.9527830433393649212092669e-01 1.2170300415936896685359381e-03
    -9.9398303667390686122473653e-01 1.3734713364273324134878784e-03
    -9.9253143365047602486583855e-01 1.5296966649272487719091185e-03
    -9.9092372353162150311334244e-01 1.6856814323170034482235469e-03
    -9.8916015915543675784959987e-01 1.8414010924821320840105709e-03
    -9.8724101788260765211191483e-01 1.9968311464057843240826884e-03
    -9.8516660154880342226135781e-01 2.1519471434930347031322384e-03
    -9.8293723641503172316902237e-01 2.3067246841547667073057948e-03
    -9.8055327311507989307415301e-01 2.4611394229790972687510475e-03
    -9.7801508659962455016767535e-01 2.6151670721914336593949546e-03
    -9.7532307607679913363796231e-01 2.7687834052614944228831728e-03
    -9.7247766494911269674616960e-01 2.9219642605862617551482074e-03
    -9.6947930074666410771300207e-01 3.0746855452115097982745962e-03
    -9.6632845505662412488590007e-01 3.2269232385712526278709333e-03
    -9.6302562344897502111251697e-01 3.3786533962332863240130010e-03
    -9.5957132539850387153990141e-01 3.5298521536436060740127285e-03
    -9.5596610420305461186529783e-01 3.6804957298652175162745337e-03
    -9.5221052689804486224289803e-01 3.8305604313083299779230106e-03
    -9.4830518416725817498758033e-01 3.9800226554498052689012866e-03
    -9.4425069024992303035048735e-01 4.1288588945403645324994102e-03
    -9.4004768284409001566359620e-01 4.2770457392982408340031952e-03
    -9.3569682300632472937707007e-01 4.4245598825883641555534176e-03
    -9.3119879504772695710812513e-01 4.5713781230861706694756919e-03
    -9.2655430642629521553743643e-01 4.7174773689252482919420650e-03
    -9.2176408763565187420851998e-01 4.8628346413281484800217314e-03
    -9.1682889209014606368697287e-01 5.0074270782196781889861192e-03
    -9.1174949600635202262566281e-01 5.1512319378220153970326933e-03
    -9.0652669828098231263879825e-01 5.2942266022310661977012813e-03
    -9.0116132036523399762018016e-01 5.4363885809734580414898097e-03
    -8.9565420613558721640856675e-01 5.5776955145435532998354766e-03
    -8.9000622176107802019373594e-01 5.7181251779199393156516429e-03
    -8.8421825556706423743236201e-01 5.8576554840608311552907495e-03
    -8.7829121789550768095722333e-01 5.9962644873777937865044763e-03
    -8.7222604096179268129418460e-01 6.1339303871872493875705423e-03
    -8.6602367870810503802658786e-01 6.2706315311392294212233800e-03
    -8.5968510665339392673445218e-01 6.4063464186228068208972530e-03
    -8.5321132173993996516969673e-01 6.5410537041476687902807896e-03
    -8.4660334217655353050702161e-01 6.6747322007013436262479189e-03
    -8.3986220727842975097843237e-01 6.8073608830814501258199556e-03
    -8.3298897730368137715117882e-01 6.9389188912025697233976196e-03
    -8.2598473328658073011609986e-01 7.0693855333770891022360239e-03
    -8.1885057686753104366772504e-01 7.1987402895696455187857232e-03
    -8.1158763011979939694384711e-01 7.3269628146244878894033370e-03
    -8.0419703537303399709657015e-01 7.4540329414654040926513368e-03
    -7.9667995503359700926182541e-01 7.5799306842676186096730007e-03
    -7.8903757140173835793461876e-01 7.7046362416011661064518812e-03
    -7.8127108648564225212851397e-01 7.8281299995453055423633160e-03
    -7.7338172181237152535970836e-01 7.9503925347734153566969795e-03
    -7.6537071823574520657018638e-01 8.0714046176079117661528173e-03
    -7.5723933574117352485188803e-01 8.1911472150447365431213953e-03
    -7.4898885324748543101947007e-01 8.3096014937468869188119669e-03
    -7.4062056840577794591951033e-01 8.4267488230065426546566698e-03
    -7.3213579739531897683235684e-01 8.5425707776753505295896929e-03
    -7.2353587471653590945663836e-01 8.6570491410623579692229512e-03
    -7.1482215298112383727868746e-01 8.7701659077991430540910400e-03
    -7.0599600269930296025933103e-01 8.8819032866717662066946559e-03
    -6.9705881206426334451009552e-01 8.9922437034189537513606538e-03
    -6.8801198673382668591358424e-01 9.1011698034962566078442947e-03
    -6.7885694960936182607014189e-01 9.2086644548055717612866644e-03
    -6.6959514061198854850687212e-01 9.3147107503897073194076839e-03
    -6.6022801645610351695125928e-01 9.4192920110915748233004763e-03
    -6.5075705042026632529683638e-01 9.5223917881774450266973986e-03
    -6.4118373211547896595163820e-01 9.6239938659241196922122796e-03
    -6.3150956725089713028609140e-01 9.7240822641692994943163342e-03
    -6.2173607739700975649554948e-01 9.8226412408250144553401029e-03
    -6.1186479974632290712577287e-01 9.9196552943535427210308697e-03
    -6.0189728687158816633484548e-01 1.0015109166205538385185925e-02
    -5.9183510648161052891680356e-01 1.0108987843219903457470110e-02
    -5.8167984117467663729428295e-01 1.0201276559985105601979782e-02
    -5.7143308818964055895150977e-01 1.0291960801161506230960718e-02
    -5.6109645915470718335882339e-01 1.0381026303764332927026537e-02
    -5.5067157983395198517229119e-01 1.0468459059407015168674526e-02
    -5.4016008987161889809414106e-01 1.0554245316504454882400310e-02
    -5.2956364253423304777612657e-01 1.0638371582435858855864019e-02
    -5.1888390445057119837457549e-01 1.0720824625666935253631706e-02
    -5.0812255534953476576731646e-01 1.0801591477830933818449211e-02
    -4.9728128779595404118651913e-01 1.0880659435768354550977399e-02
    -4.8636180692438263362120665e-01 1.0958016063524910269078028e-02
    -4.7536583017090749958555307e-01 1.1033649194307500185363580e-02
    -4.6429508700303090407146556e-01 1.1107546932397786842994236e-02
    -4.5315131864765362257330139e-01 1.1179697655023209876268275e-02
    -4.4193627781721228631184317e-01 1.1250090014185035874882956e-02
    -4.3065172843401033908605768e-01 1.1318712938443161591939301e-02
    -4.1929944535278468320527168e-01 1.1385555634657479676108416e-02
    -4.0788121408155358915692545e-01 1.1450607589685442636029400e-02
    -3.9639883050078744686217647e-01 1.1513858572035563093693966e-02
    -3.8485410058095026464286548e-01 1.1575298633476673385023226e-02
    -3.7324884009845393784132739e-01 1.1634918110602590624047536e-02
    -3.6158487435006858579100708e-01 1.1692707626351968655531444e-02
    -3.4986403786583664121678794e-01 1.1748658091483194648718680e-02
    -3.3808817412053493445256436e-01 1.1802760706003905324945613e-02
    -3.2625913524372973650855556e-01 1.1855006960555096798271002e-02
    -3.1437878172846911439819451e-01 1.1905388637749496696938145e-02
    -3.0244898213866316938336354e-01 1.1953897813463979990511454e-02
    -2.9047161281519084941038500e-01 1.2000526858085929643449319e-02
    -2.7844855758078806973188080e-01 1.2045268437713196885141542e-02
    -2.6638170744375311294049880e-01 1.2088115515307610789430548e-02
    -2.5427296030052870534632348e-01 1.2129061351801789098159290e-02
    -2.4212422063719571396767094e-01 1.2168099507159043856652225e-02
    -2.2993739922993181035160148e-01 1.2205223841386302063849456e-02
    -2.1771441284448234121384758e-01 1.2240428515499800210530879e-02
    -2.0545718393468664908496635e-01 1.2273707992443452091668732e-02
    -1.9316764034011163486681539e-01 1.2305057037959768462265231e-02
    -1.8084771498283738755397110e-01 1.2334470721413037103131316e-02
    -1.6849934556344323133281193e-01 1.2361944416564876053632460e-02
    -1.5612447425624337293825761e-01 1.2387473802301837352835001e-02
    -1.4372504740381872312404710e-01 1.2411054863315050880712143e-02
    -1.3130301521089041139056519e-01 1.2432683890731746970126359e-02
    -1.1886033143758943653178761e-01 1.2452357482698607482607400e-02
    -1.0639895309216537699903427e-01 1.2470072544916794426983486e-02
    -9.3920840123185161951724353e-02 1.2485826291128645712524836e-02
    -8.1427955111268157661896794e-02 1.2499616243555916605956213e-02
    -6.8922262960407404408513798e-02 1.2511440233289441992248214e-02
    -5.6405730588925909185782359e-02 1.2521296400630321787872390e-02
    -4.3880326620117995894965190e-02 1.2529183195382365828551841e-02
    -3.1348021072617915372404696e-02 1.2535099377095959560790561e-02
    -1.8810785050553550934449021e-02 1.2539044015263126757853129e-02
    -6.2705904335270835209259488e-03 1.2541016489463895078326772e-02
    6.2705904335270835209259488e-03 1.2541016489463895078326772e-02
    1.8810785050553550934449021e-02 1.2539044015263126757853129e-02
    3.1348021072617915372404696e-02 1.2535099377095959560790561e-02
    4.3880326620117995894965190e-02 1.2529183195382365828551841e-02
    5.6405730588925909185782359e-02 1.2521296400630321787872390e-02
    6.8922262960407404408513798e-02 1.2511440233289441992248214e-02
    8.1427955111268157661896794e-02 1.2499616243555916605956213e-02
    9.3920840123185161951724353e-02 1.2485826291128645712524836e-02
    1.0639895309216537699903427e-01 1.2470072544916794426983486e-02
    1.1886033143758943653178761e-01 1.2452357482698607482607400e-02
    1.3130301521089041139056519e-01 1.2432683890731746970126359e-02
    1.4372504740381872312404710e-01 1.2411054863315050880712143e-02
    1.5612447425624337293825761e-01 1.2387473802301837352835001e-02
    1.6849934556344323133281193e-01 1.2361944416564876053632460e-02
    1.8084771498283738755397110e-01 1.2334470721413037103131316e-02
    1.9316764034011163486681539e-01 1.2305057037959768462265231e-02
    2.0545718393468664908496635e-01 1.2273707992443452091668732e-02
    2.1771441284448234121384758e-01 1.2240428515499800210530879e-02
    2.2993739922993181035160148e-01 1.2205223841386302063849456e-02
    2.4212422063719571396767094e-01 1.2168099507159043856652225e-02
    2.5427296030052870534632348e-01 1.2129061351801789098159290e-02
    2.6638170744375311294049880e-01 1.2088115515307610789430548e-02
    2.7844855758078806973188080e-01 1.2045268437713196885141542e-02
    2.9047161281519084941038500e-01 1.2000526858085929643449319e-02
    3.0244898213866316938336354e-01 1.1953897813463979990511454e-02
    3.1437878172846911439819451e-01 1.1905388637749496696938145e-02
    3.2625913524372973650855556e-01 1.1855006960555096798271002e-02
    3.3808817412053493445256436e-01 1.1802760706003905324945613e-02
    3.4986403786583664121678794e-01 1.1748658091483194648718680e-02
    3.6158487435006858579100708e-01 1.1692707626351968655531444e-02
    3.7324884009845393784132739e-01 1.1634918110602590624047536e-02
    3.8485410058095026464286548e-01 1.1575298633476673385023226e-02
    3.9639883050078744686217647e-01 1.1513858572035563093693966e-02
    4.0788121408155358915692545e-01 1.1450607589685442636029400e-02
    4.1929944535278468320527168e-01 1.1385555634657479676108416e-02
    4.3065172843401033908605768e-01 1.1318712938443161591939301e-02
    4.4193627781721228631184317e-01 1.1250090014185035874882956e-02
    4.5315131864765362257330139e-01 1.1179697655023209876268275e-02
    4.6429508700303090407146556e-01 1.1107546932397786842994236e-02
    4.7536583017090749958555307e-01 1.1033649194307500185363580e-02
    4.8636180692438263362120665e-01 1.0958016063524910269078028e-02
    4.9728128779595404118651913e-01 1.0880659435768354550977399e-02
    5.0812255534953476576731646e-01 1.0801591477830933818449211e-02
    5.1888390445057119837457549e-01 1.0720824625666935253631706e-02
    5.2956364253423304777612657e-01 1.0638371582435858855864019e-02
    5.4016008987161889809414106e-01 1.0554245316504454882400310e-02
    5.5067157983395198517229119e-01 1.0468459059407015168674526e-02
    5.6109645915470718335882339e-01 1.0381026303764332927026537e-02
    5.7143308818964055895150977e-01 1.0291960801161506230960718e-02
    5.8167984117467663729428295e-01 1.0201276559985105601979782e-02
    5.9183510648161052891680356e-01 1.0108987843219903457470110e-02
    6.0189728687158816633484548e-01 1.0015109166205538385185925e-02
    6.1186479974632290712577287e-01 9.9196552943535427210308697e-03
    6.2173607739700975649554948e-01 9.8226412408250144553401029e-03
    6.3150956725089713028609140e-01 9.7240822641692994943163342e-03
    6.4118373211547896595163820e-01 9.6239938659241196922122796e-03
    6.5075705042026632529683638e-01 9.5223917881774450266973986e-03
    6.6022801645610351695125928e-01 9.4192920110915748233004763e-03
    6.6959514061198854850687212e-01 9.3147107503897073194076839e-03
    6.7885694960936182607014189e-01 9.2086644548055717612866644e-03
    6.8801198673382668591358424e-01 9.1011698034962566078442947e-03
    6.9705881206426334451009552e-01 8.9922437034189537513606538e-03
    7.0599600269930296025933103e-01 8.8819032866717662066946559e-03
    7.1482215298112383727868746e-01 8.7701659077991430540910400e-03
    7.2353587471653590945663836e-01 8.6570491410623579692229512e-03
    7.3213579739531897683235684e-01 8.5425707776753505295896929e-03
    7.4062056840577794591951033e-01 8.4267488230065426546566698e-03
    7.4898885324748543101947007e-01 8.3096014937468869188119669e-03
    7.5723933574117352485188803e-01 8.1911472150447365431213953e-03
    7.6537071823574520657018638e-01 8.0714046176079117661528173e-03
    7.7338172181237152535970836e-01 7.9503925347734153566969795e-03
    7.8127108648564225212851397e-01 7.8281299995453055423633160e-03
    7.8903757140173835793461876e-01 7.7046362416011661064518812e-03
    7.9667995503359700926182541e-01 7.5799306842676186096730007e-03
    8.0419703537303399709657015e-01 7.4540329414654040926513368e-03
    8.1158763011979939694384711e-01 7.3269628146244878894033370e-03
    8.1885057686753104366772504e-01 7.1987402895696455187857232e-03
    8.2598473328658073011609986e-01 7.0693855333770891022360239e-03
    8.3298897730368137715117882e-01 6.9389188912025697233976196e-03
    8.3986220727842975097843237e-01 6.8073608830814501258199556e-03
    8.4660334217655353050702161e-01 6.6747322007013436262479189e-03
    8.5321132173993996516969673e-01 6.5410537041476687902807896e-03
    8.5968510665339392673445218e-01 6.4063464186228068208972530e-03
    8.6602367870810503802658786e-01 6.2706315311392294212233800e-03
    8.7222604096179268129418460e-01 6.1339303871872493875705423e-03
    8.7829121789550768095722333e-01 5.9962644873777937865044763e-03
    8.8421825556706423743236201e-01 5.8576554840608311552907495e-03
    8.9000622176107802019373594e-01 5.7181251779199393156516429e-03
    8.9565420613558721640856675e-01 5.5776955145435532998354766e-03
    9.0116132036523399762018016e-01 5.4363885809734580414898097e-03
    9.0652669828098231263879825e-01 5.2942266022310661977012813e-03
    9.1174949600635202262566281e-01 5.1512319378220153970326933e-03
    9.1682889209014606368697287e-01 5.0074270782196781889861192e-03
    9.2176408763565187420851998e-01 4.8628346413281484800217314e-03
    9.2655430642629521553743643e-01 4.7174773689252482919420650e-03
    9.3119879504772695710812513e-01 4.5713781230861706694756919e-03
    9.3569682300632472937707007e-01 4.4245598825883641555534176e-03
    9.4004768284409001566359620e-01 4.2770457392982408340031952e-03
    9.4425069024992303035048735e-01 4.1288588945403645324994102e-03
    9.4830518416725817498758033e-01 3.9800226554498052689012866e-03
    9.5221052689804486224289803e-01 3.8305604313083299779230106e-03
    9.5596610420305461186529783e-01 3.6804957298652175162745337e-03
    9.5957132539850387153990141e-01 3.5298521536436060740127285e-03
    9.6302562344897502111251697e-01 3.3786533962332863240130010e-03
    9.6632845505662412488590007e-01 3.2269232385712526278709333e-03
    9.6947930074666410771300207e-01 3.0746855452115097982745962e-03
    9.7247766494911269674616960e-01 2.9219642605862617551482074e-03
    9.7532307607679913363796231e-01 2.7687834052614944228831728e-03
    9.7801508659962455016767535e-01 2.6151670721914336593949546e-03
    9.8055327311507989307415301e-01 2.4611394229790972687510475e-03
    9.8293723641503172316902237e-01 2.3067246841547667073057948e-03
    9.8516660154880342226135781e-01 2.1519471434930347031322384e-03
    9.8724101788260765211191483e-01 1.9968311464057843240826884e-03
    9.8916015915543675784959987e-01 1.8414010924821320840105709e-03
    9.9092372353162150311334244e-01 1.6856814323170034482235469e-03
    9.9253143365047602486583855e-01 1.5296966649272487719091185e-03
    9.9398303667390686122473653e-01 1.3734713364273324134878784e-03
    9.9527830433393649212092669e-01 1.2170300415936896685359381e-03
    9.9641703298482942052771705e-01 1.0603974326420877692667144e-03
    9.9739904367225573622590673e-01 9.0359824806269411524850543e-04
    9.9822418225691300630586511e-01 7.4665740556624201132440710e-04
    9.9889231972585434959910344e-01 5.8960034219324988433164059e-04
    9.9940335329301610567398484e-01 4.3245460060038212139060798e-04
    9.9975721221158908580406433e-01 2.7526103865831749294840192e-04
    9.9995391943564204684236074e-01 1.1825670068171721075345887e-04];

X=xw(:,1); W=xw(:,2);







function [X,W]=gauss_legendre_700

xw=[ -9.9999410721789827594108147e-01 1.5122766974076512616403586e-05
    -9.9996895141230646153474027e-01 3.5202627719042746520707737e-05
    -9.9992369471822484250367324e-01 5.5311513823237843051703472e-05
    -9.9985832799377194479717446e-01 7.5421869998298820669119236e-05
    -9.9977285132240067966336028e-01 9.5531227691402493119691985e-05
    -9.9966726612928136219693442e-01 1.1563881924509692518312359e-04
    -9.9954157444010394151234777e-01 1.3574414394147891441307074e-04
    -9.9939577874291152248531489e-01 1.5584676466203183911249375e-04
    -9.9922988195134376798733911e-01 1.7594626390712142906605497e-04
    -9.9904388739237171002116611e-01 1.9604223142283443213733374e-04
    -9.9883779880148382268600926e-01 2.1613426002684216816612794e-04
    -9.9861162032053651937957284e-01 2.3622194399989265996796239e-04
    -9.9836535649667978997712225e-01 2.5630487840044917065632224e-04
    -9.9809901228175612608595202e-01 2.7638265874956383354871248e-04
    -9.9781259303192360032852548e-01 2.9645488087773699554389184e-04
    -9.9750610450738885770505249e-01 3.1652114084790224394702629e-04
    -9.9717955287219695037492784e-01 3.3658103491638908383309925e-04
    -9.9683294469404826187997060e-01 3.5663415951382569256053467e-04
    -9.9646628694412919813316876e-01 3.7668011123697140796542926e-04
    -9.9607958699694720827721994e-01 3.9671848684678294536670728e-04
    -9.9567285263016636065458442e-01 4.1674888327014565159514548e-04
    -9.9524609202443881095234701e-01 4.3677089760383425815978842e-04
    -9.9479931376323316172261002e-01 4.5678412711985843163917109e-04
    -9.9433252683265593852013353e-01 4.7678816927169146095213947e-04
    -9.9384574062126618265722300e-01 4.9678262170107492703280405e-04
    -9.9333896491988438182119125e-01 5.1676708224520316852651503e-04
    -9.9281220992139318504143830e-01 5.3674114894416698022527878e-04
    -9.9226548622053001302845132e-01 5.5670442004856959122838234e-04
    -9.9169880481367378433077420e-01 5.7665649402726456163731372e-04
    -9.9111217709862164948475538e-01 5.9659696957518036121437266e-04
    -9.9050561487436006302687019e-01 6.1652544562120089287338143e-04
    -9.8987913034082575247651903e-01 6.3644152133608743269338470e-04
    -9.8923273609866091415909750e-01 6.5634479614043493914121719e-04
    -9.8856644514895852804414744e-01 6.7623486971263996325731682e-04
    -9.8788027089300134431226752e-01 6.9611134199689035130520498e-04
    -9.8717422713199121098170963e-01 7.1597381321115885054878758e-04
    -9.8644832806677229530833984e-01 7.3582188385520343708312563e-04
    -9.8570258829754442420068017e-01 7.5565515471857100469194046e-04
    -9.8493702282356965227450019e-01 7.7547322688860127054555349e-04
    -9.8415164704287072527932878e-01 7.9527570175842797039339471e-04
    -9.8334647675192077276307145e-01 8.1506218103498189690020448e-04
    -9.8252152814532589530926998e-01 8.3483226674698667382767958e-04
    -9.8167681781549931407937493e-01 8.5458556125295192813090539e-04
    -9.8081236275232719368233347e-01 8.7432166724916331786848778e-04
    -9.7992818034282747063912211e-01 8.9404018777766638016030187e-04
    -9.7902428837079980006308233e-01 9.1374072623424116342677470e-04
    -9.7810070501646784180138638e-01 9.3342288637637761856968854e-04
    -9.7715744885611366399302824e-01 9.5308627233123993129082496e-04
    -9.7619453886170381995412981e-01 9.7273048860362402393692216e-04
    -9.7521199440050809759128470e-01 9.9235514008391115421681139e-04
    -9.7420983523470983111991472e-01 1.0119598320560106718968285e-03
    -9.7318808152100833019915171e-01 1.0315441702052981134268839e-03
    -9.7214675381021331546094189e-01 1.0511077606265416956032865e-03
    -9.7108587304683191554488531e-01 1.0706502098318268576676582e-03
    -9.7000546056864733746039064e-01 1.0901711247584676845595597e-03
    -9.6890553810628932129844770e-01 1.1096701127769103070730417e-03
    -9.6778612778279737849373987e-01 1.1291467816986319652045045e-03
    -9.6664725211317581443637437e-01 1.1486007397840248926701445e-03
    -9.6548893400394086850013764e-01 1.1680315957502765151360125e-03
    -9.6431119675265963842036854e-01 1.1874389587792388064263482e-03
    -9.6311406404748201026677634e-01 1.2068224385252889711145352e-03
    -9.6189755996666392867666673e-01 1.2261816451231773335284192e-03
    -9.6066170897808289552699534e-01 1.2455161891958719741491102e-03
    -9.5940653593874647420136625e-01 1.2648256818623868861367621e-03
    -9.5813206609429157900592600e-01 1.2841097347456062287285317e-03
    -9.5683832507847732529171481e-01 1.3033679599800925344443847e-03
    -9.5552533891266877574821592e-01 1.3225999702198892784416051e-03
    -9.5419313400531458047026945e-01 1.3418053786463130395284482e-03
    -9.5284173715141451399546213e-01 1.3609837989757285307834689e-03
    -9.5147117553198179429330139e-01 1.3801348454673253648983255e-03
    -9.5008147671349552076947020e-01 1.3992581329308674976258375e-03
    -9.4867266864734633990963175e-01 1.4183532767344489596006429e-03
    -9.4724477966927489447357402e-01 1.4374198928122235673077167e-03
    -9.4579783849880105783825002e-01 1.4564575976721339835295854e-03
    -9.4433187423864661802497267e-01 1.4754660084036210453140026e-03
    -9.4284691637415030118773984e-01 1.4944447426853291888138031e-03
    -9.4134299477267413536196727e-01 1.5133934187927951774133017e-03
    -9.3982013968300293083046881e-01 1.5323116556061207836625382e-03
    -9.3827838173473598892826431e-01 1.5511990726176452543710882e-03
    -9.3671775193767115030851755e-01 1.5700552899395952242966867e-03
    -9.3513828168118051653578959e-01 1.5888799283117253056951679e-03
    -9.3354000273357984340805160e-01 1.6076726091089474021678107e-03
    -9.3192294724148894147219835e-01 1.6264329543489433257952292e-03
    -9.3028714772918563724601881e-01 1.6451605866997726437817029e-03
    -9.2863263709795118572287720e-01 1.6638551294874545042473679e-03
    -9.2695944862540868847133879e-01 1.6825162067035522809460568e-03
    -9.2526761596485385119592593e-01 1.7011434430127274266414394e-03
    -9.2355717314457808075900402e-01 1.7197364637602950612066399e-03
    -9.2182815456718381064149526e-01 1.7382948949797596104038799e-03
    -9.2008059500889272097623461e-01 1.7568183634003315795463207e-03
    -9.1831452961884640906475852e-01 1.7753064964544380387873046e-03
    -9.1652999391839951037752598e-01 1.7937589222852166285365749e-03
    -9.1472702380040482594836249e-01 1.8121752697539956870886879e-03
    -9.1290565552849134434154621e-01 1.8305551684477561636549270e-03
    -9.1106592573633582432535150e-01 1.8488982486865846409374026e-03
    -9.0920787142692471860527803e-01 1.8672041415311048905001368e-03
    -9.0733152997181054644215692e-01 1.8854724787899042239697200e-03
    -9.0543693911036049470908438e-01 1.9037028930269346907455663e-03
    -9.0352413694899691432027566e-01 1.9218950175689025663267051e-03
    -9.0159316196043082225486387e-01 1.9400484865126482996594559e-03
    -8.9964405298288829815334111e-01 1.9581629347325008733138318e-03
    -8.9767684921932955344203720e-01 1.9762379978876267426490809e-03
    -8.9569159023665978480721606e-01 1.9942733124293514530844806e-03
    -8.9368831596493392144253676e-01 2.0122685156084758363592702e-03
    -8.9166706669655348971303965e-01 2.0302232454825746427173883e-03
    -8.8962788308545570625796017e-01 2.0481371409232698027613750e-03
    -8.8757080614629635384460471e-01 2.0660098416235002198593218e-03
    -8.8549587725362433054954181e-01 2.0838409881047633732953361e-03
    -8.8340313814105009271315794e-01 2.1016302217243491151643653e-03
    -8.8129263090040554917692361e-01 2.1193771846825552526705216e-03
    -8.7916439798089818236093151e-01 2.1370815200298736401263167e-03
    -8.7701848218825650960184248e-01 2.1547428716741862456118817e-03
    -8.7485492668386966030880103e-01 2.1723608843879096096107784e-03
    -8.7267377498391918155817848e-01 2.1899352038151527477527480e-03
    -8.7047507095850373826095847e-01 2.2074654764788399254060725e-03
    -8.6825885883075637483585751e-01 2.2249513497878173860777817e-03
    -8.6602518317595578167811254e-01 2.2423924720439488041112686e-03
    -8.6377408892062879086637395e-01 2.2597884924491899207021905e-03
    -8.6150562134164809791059270e-01 2.2771390611126388937857090e-03
    -8.5921982606531999149268586e-01 2.2944438290575784079872168e-03
    -8.5691674906646797538201099e-01 2.3117024482284974015722234e-03
    -8.5459643666750717549973615e-01 2.3289145714980807010119346e-03
    -8.5225893553751330689038923e-01 2.3460798526742051607618667e-03
    -8.4990429269128375810993248e-01 2.3631979465068889655066098e-03
    -8.4753255548839157018647938e-01 2.3802685086952457028941499e-03
    -8.4514377163223375344358601e-01 2.3972911958944072112476231e-03
    -8.4273798916907072253934530e-01 2.4142656657224238758718826e-03
    -8.4031525648706084652417303e-01 2.4311915767671579864661613e-03
    -8.3787562231528678324821158e-01 2.4480685885931484715993811e-03
    -8.3541913572277581856440065e-01 2.4648963617484556838654886e-03
    -8.3294584611751176783656092e-01 2.4816745577714927409318513e-03
    -8.3045580324544276962228651e-01 2.4984028391978295446929081e-03
    -8.2794905718947975348243062e-01 2.5150808695669811879125355e-03
    -8.2542565836848913463086319e-01 2.5317083134291798809933915e-03
    -8.2288565753627918031298805e-01 2.5482848363521134853193484e-03
    -8.2032910578057882666769274e-01 2.5648101049276579750657579e-03
    -8.1775605452201027834036040e-01 2.5812837867785843118384470e-03
    -8.1516655551305472471312896e-01 2.5977055505652392984605736e-03
    -8.1256066083701139479700259e-01 2.6140750659922112539290673e-03
    -8.0993842290695006180811788e-01 2.6303920038149761727319653e-03
    -8.0729989446465699742816469e-01 2.6466560358465143938666575e-03
    -8.0464512857957359859284452e-01 2.6628668349639194636024886e-03
    -8.0197417864772979623211313e-01 2.6790240751149692680077585e-03
    -7.9928709839066958142694830e-01 2.6951274313246884918593427e-03
    -7.9658394185437009227257477e-01 2.7111765797018880924662554e-03
    -7.9386476340815570473807838e-01 2.7271711974456679106193491e-03
    -7.9112961774360357480873063e-01 2.7431109628519253530731792e-03
    -7.8837855987344362951318999e-01 2.7589955553198146354088038e-03
    -7.8561164513045289581327779e-01 2.7748246553581956165557276e-03
    -7.8282892916634128077646437e-01 2.7905979445920674544834306e-03
    -7.8003046795063346596776910e-01 2.8063151057689645316572946e-03
    -7.7721631776954280823588306e-01 2.8219758227653406711110673e-03
    -7.7438653522483846813884156e-01 2.8375797805929251632628407e-03
    -7.7154117723270820849990059e-01 2.8531266654050636139006158e-03
    -7.6868030102261264424612364e-01 2.8686161645030232303366091e-03
    -7.6580396413613460726566018e-01 2.8840479663422885665824413e-03
    -7.6291222442582184992687644e-01 2.8994217605388225740548069e-03
    -7.6000514005402330930394328e-01 2.9147372378753129071315975e-03
    -7.5708276949171915415348622e-01 2.9299940903073887384089335e-03
    -7.5414517151734494770920492e-01 2.9451920109698193593616988e-03
    -7.5119240521560926016064741e-01 2.9603306941826876275136726e-03
    -7.4822452997630561899455870e-01 2.9754098354575278517764048e-03
    -7.4524160549311724288656933e-01 2.9904291315034658461924355e-03
    -7.4224369176241744572308789e-01 3.0053882802333021714846772e-03
    -7.3923084908206193599511380e-01 3.0202869807696049175849140e-03
    -7.3620313805017678632225397e-01 3.0351249334507495770962837e-03
    -7.3316061956393907550477707e-01 3.0499018398369498114575116e-03
    -7.3010335481835231252745189e-01 3.0646174027162604615315278e-03
    -7.2703140530501531024754058e-01 3.0792713261105541036610767e-03
    -7.2394483281088617410148345e-01 3.0938633152814754880000070e-03
    -7.2084369941703840822810889e-01 3.1083930767363695223115538e-03
    -7.1772806749741302478895477e-01 3.1228603182341862706805635e-03
    -7.1459799971756421399504688e-01 3.1372647487913521251179372e-03
    -7.1145355903339813075092479e-01 3.1516060786876305688242184e-03
    -7.0829480868990712938426668e-01 3.1658840194719460765793162e-03
    -7.0512181221989789214887878e-01 3.1800982839681828616418269e-03
    -7.0193463344271289638953704e-01 3.1942485862809654080518573e-03
    -6.9873333646294777388163766e-01 3.2083346418014073442304124e-03
    -6.9551798566916078758737285e-01 3.2223561672128338620457733e-03
    -6.9228864573257897774283265e-01 3.2363128804964850539216403e-03
    -6.8904538160579675842853931e-01 3.2502045009371819533905423e-03
    -6.8578825852146896302485857e-01 3.2640307491289782641785688e-03
    -6.8251734199100011490912721e-01 3.2777913469807713232884172e-03
    -6.7923269780322581556930572e-01 3.2914860177219108956780946e-03
    -6.7593439202308958080323009e-01 3.3051144859077511567457819e-03
    -6.7262249099031456989195021e-01 3.3186764774251957359207932e-03
    -6.6929706131806887547952556e-01 3.3321717194982150046789204e-03
    -6.6595816989162581744920999e-01 3.3455999406933256343221750e-03
    -6.6260588386701835261760607e-01 3.3589608709250575770133285e-03
    -6.5924027066968915455902334e-01 3.3722542414613772450426854e-03
    -6.5586139799313358800247897e-01 3.3854797849291002817539997e-03
    -6.5246933379753857540350737e-01 3.3986372353192679032773960e-03
    -6.4906414630841602342314900e-01 3.4117263279924928826014607e-03
    -6.4564590401522947704648914e-01 3.4247467996842899211340505e-03
    -6.4221467567001810916593740e-01 3.4376983885103609174527239e-03
    -6.3877053028601249451412514e-01 3.4505808339718719961186633e-03
    -6.3531353713624749701693872e-01 3.4633938769606828315950064e-03
    -6.3184376575216827376380024e-01 3.4761372597645729363990963e-03
    -6.2836128592223272626426933e-01 3.4888107260724063665713324e-03
    -6.2486616769050606912117019e-01 3.5014140209793051007614029e-03
    -6.2135848135525451052529888e-01 3.5139468909917625713545952e-03
    -6.1783829746752882972060661e-01 3.5264090840327480882998579e-03
    -6.1430568682974717731326564e-01 3.5388003494467757011066045e-03
    -6.1076072049426977095265556e-01 3.5511204380049496420745925e-03
    -6.0720346976197092647709042e-01 3.5633691019099694372029852e-03
    -6.0363400618080398363218819e-01 3.5755460948011198382689990e-03
    -6.0005240154436256805325911e-01 3.5876511717592247593944244e-03
    -5.9645872789043674622178060e-01 3.5996840893115725906747926e-03
    -5.9285305749956462850747130e-01 3.6116446054368115878285828e-03
    -5.8923546289357675576070505e-01 3.6235324795698235778029783e-03
    -5.8560601683413981977111007e-01 3.6353474726065564646970696e-03
    -5.8196479232129150194197109e-01 3.6470893469088324495563924e-03
    -5.7831186259197264742937250e-01 3.6587578663091358671666242e-03
    -5.7464730111855444327773057e-01 3.6703527961153554863560355e-03
    -5.7097118160736048952941246e-01 3.6818739031155077283397059e-03
    -5.6728357799718420739765179e-01 3.6933209555824290937220411e-03
    -5.6358456445780258370348292e-01 3.7046937232784317266254259e-03
    -5.5987421538848236579610784e-01 3.7159919774599385958180342e-03
    -5.5615260541648581238405313e-01 3.7272154908820887518616782e-03
    -5.5241980939556856178285216e-01 3.7383640378032957467258068e-03
    -5.4867590240447439153825826e-01 3.7494373939897947776989096e-03
    -5.4492095974542453795663732e-01 3.7604353367201568715538773e-03
    -5.4115505694260446212240367e-01 3.7713576447897614353499396e-03
    -5.3737826974064339946579594e-01 3.7822040985152531947233179e-03
    -5.3359067410309157786230116e-01 3.7929744797389488251970757e-03
    -5.2979234621089221768386324e-01 3.8036685718332353435544579e-03
    -5.2598336246084964606950507e-01 3.8142861597049272995296754e-03
    -5.2216379946409208212543263e-01 3.8248270297995823341352661e-03
    -5.1833373404453042532225027e-01 3.8352909701058102327764221e-03
    -5.1449324323731515651303425e-01 3.8456777701595255998523459e-03
    -5.1064240428728369103339446e-01 3.8559872210481914260593594e-03
    -5.0678129464740995224758535e-01 3.8662191154150088792662565e-03
    -5.0290999197724384206509285e-01 3.8763732474631023249000972e-03
    -4.9902857414135165514679215e-01 3.8864494129596405952442595e-03
    -4.9513711920774744479345486e-01 3.8964474092399651976303954e-03
    -4.9123570544632555456843193e-01 3.9063670352116513020956923e-03
    -4.8732441132728337995771994e-01 3.9162080913585578870184101e-03
    -4.8340331551954623945377421e-01 3.9259703797448440576456363e-03
    -4.7947249688918214260979767e-01 3.9356537040189367923637143e-03
    -4.7553203449781783485050823e-01 3.9452578694174995563304087e-03
    -4.7158200760104723725518738e-01 3.9547826827693380313810678e-03
    -4.6762249564683749936122581e-01 3.9642279524992815598061213e-03
    -4.6365357827393199885435138e-01 3.9735934886320481082555567e-03
    -4.5967533531024740156567532e-01 3.9828791027960563225773782e-03
    -4.5568784677126900062305026e-01 3.9920846082272185006978305e-03
    -4.5169119285843972733118790e-01 4.0012098197726945342234117e-03
    -4.4768545395754927307407911e-01 4.0102545538946180944672726e-03
    -4.4367071063711610579005651e-01 4.0192186286737872566443563e-03
    -4.3964704364676626680008553e-01 4.0281018638133212969587582e-03
    -4.3561453391561111292418218e-01 4.0369040806422940709241587e-03
    -4.3157326255061706499205343e-01 4.0456251021193162173417157e-03
    -4.2752331083497724373287951e-01 4.0542647528361086886605769e-03
    -4.2346476022647383530284060e-01 4.0628228590210199028254223e-03
    -4.1939769235584140050221436e-01 4.0712992485425281499744621e-03
    -4.1532218902512429981044306e-01 4.0796937509127023657740096e-03
    -4.1123833220603150939709280e-01 4.0880061972906264755600603e-03
    -4.0714620403828893913100728e-01 4.0962364204858055238833714e-03
    -4.0304588682798497822901140e-01 4.1043842549615102213711459e-03
    -3.9893746304591748419454689e-01 4.1124495368381310325678335e-03
    -3.9482101532593349979549657e-01 4.1204321038964550685812505e-03
    -3.9069662646326763777437918e-01 4.1283317955809334387495824e-03
    -3.8656437941287757897868005e-01 4.1361484530029147752006047e-03
    -3.8242435728777413039836119e-01 4.1438819189438379914092181e-03
    -3.7827664335735089462531278e-01 4.1515320378583972851793682e-03
    -3.7412132104570755553041295e-01 4.1590986558776646409008038e-03
    -3.6995847392997532887548573e-01 4.1665816208121863109536953e-03
    -3.6578818573863247642918850e-01 4.1739807821550541436228521e-03
    -3.6161054034982176297319256e-01 4.1812959910849161956902797e-03
    -3.5742562178966352792741645e-01 4.1885271004689873450277382e-03
    -3.5323351423056509323927799e-01 4.1956739648659861774415702e-03
    -3.4903430198952806184919950e-01 4.2027364405290840165818622e-03
    -3.4482806952645367326582004e-01 4.2097143854087862996360414e-03
    -3.4061490144244083166924497e-01 4.2166076591557818606381680e-03
    -3.3639488247808746468336949e-01 4.2234161231237809380756332e-03
    -3.3216809751178510978775194e-01 4.2301396403722889977272459e-03
    -3.2793463155801227948415999e-01 4.2367780756693727492456780e-03
    -3.2369456976562277494835485e-01 4.2433312954943732536738921e-03
    -3.1944799741613549848295861e-01 4.2497991680405956121946431e-03
    -3.1519499992201699400951043e-01 4.2561815632179518173461297e-03
    -3.1093566282496531982815213e-01 4.2624783526555984000672161e-03
    -3.0667007179418792617298095e-01 4.2686894097044812690366911e-03
    -3.0239831262468003236776326e-01 4.2748146094399265201846383e-03
    -2.9812047123549778593343262e-01 4.2808538286641176218161320e-03
    -2.9383663366803114413983167e-01 4.2868069459085916816931672e-03
    -2.8954688608427203710959930e-01 4.2926738414366845397740491e-03
    -2.8525131476508286398896530e-01 4.2984543972459246866102411e-03
    -2.8095000610846093680450508e-01 4.3041484970704167734023571e-03
    -2.7664304662779876098355203e-01 4.3097560263831791518840397e-03
    -2.7233052295014642529835669e-01 4.3152768723984493218215341e-03
    -2.6801252181446771905015680e-01 4.3207109240739416736176715e-03
    -2.6368913006989669334245718e-01 4.3260580721130982920219488e-03
    -2.5936043467398955941760619e-01 4.3313182089672729729867839e-03
    -2.5502652269097653148222093e-01 4.3364912288378961585655347e-03
    -2.5068748129000989477432881e-01 4.3415770276786129835966399e-03
    -2.4634339774341257323087007e-01 4.3465755031973658112365300e-03
    -2.4199435942492036888396001e-01 4.3514865548584516150021351e-03
    -2.3764045380792672701453228e-01 4.3563100838845524725995162e-03
    -2.3328176846372164487952716e-01 4.3610459932587122833247406e-03
    -2.2891839105973099677271421e-01 4.3656941877262944035065217e-03
    -2.2455040935775436028443153e-01 4.3702545737969002506706495e-03
    -2.2017791121219731920177765e-01 4.3747270597462402028088313e-03
    -2.1580098456830743014478458e-01 4.3791115556179828136040832e-03
    -2.1141971746040272295275031e-01 4.3834079732255719352718160e-03
    -2.0703419801010256029449863e-01 4.3876162261539865955262130e-03
    -2.0264451442455458374247712e-01 4.3917362297614904662057533e-03
    -1.9825075499466160433570394e-01 4.3957679011813292901944550e-03
    -1.9385300809330410776176734e-01 4.3997111593233918791501225e-03
    -1.8945136217356448038451333e-01 4.4035659248758589681682629e-03
    -1.8504590576694612824582009e-01 4.4073321203067757426130591e-03
    -1.8063672748159484426899724e-01 4.4110096698656278343952941e-03
    -1.7622391600051459659148634e-01 4.4145984995848618070990454e-03
    -1.7180756007978578714379125e-01 4.4180985372813657424684308e-03
    -1.6738774854677854198037323e-01 4.4215097125579238060422149e-03
    -1.6296457029836675101108767e-01 4.4248319568046439245745383e-03
    -1.5853811429913960973081544e-01 4.4280652032003186766018210e-03
    -1.5410846957961307768236736e-01 4.4312093867137835809244528e-03
    -1.4967572523443883891758333e-01 4.4342644441052172718520374e-03
    -1.4523997042061242979116287e-01 4.4372303139274199904051876e-03
    -1.4080129435567845241905616e-01 4.4401069365270400338130408e-03
    -1.3635978631593861920556776e-01 4.4428942540457871945847046e-03
    -1.3191553563465446829994221e-01 4.4455922104215846168973059e-03
    -1.2746863170025130029827665e-01 4.4482007513896954994936372e-03
    -1.2301916395451978347708177e-01 4.4507198244838307166215685e-03
    -1.1856722189081722595993540e-01 4.4531493790371792437787768e-03
    -1.1411289505226755724986987e-01 4.4554893661834567980539745e-03
    -1.0965627302996099057263990e-01 4.4577397388578582013152207e-03
    -1.0519744546115022965082630e-01 4.4599004517980279579947300e-03
    -1.0073650202744940940213780e-01 4.4619714615449545050407387e-03
    -9.6273532453029789235365854e-02 4.4639527264438557882519909e-03
    -9.1808626502814461023405102e-02 4.4658442066450127969079453e-03
    -8.7341873980673945254693535e-02 4.4676458641045675365677248e-03
    -8.2873364727618248348761654e-02 4.4693576625852889094847598e-03
    -7.8403188619992875141306854e-02 4.4709795676573134415310307e-03
    -7.3931435567667125319246679e-02 4.4725115466988200896292227e-03
    -6.9458195512229342627463780e-02 4.4739535688966955082057630e-03
    -6.4983558425176432549186245e-02 4.4753056052471568149186965e-03
    -6.0507614306102325274494547e-02 4.4765676285563240494047577e-03
    -5.6030453180892419995373643e-02 4.4777396134407830910473258e-03
    -5.1552165099906635536974875e-02 4.4788215363280817898905539e-03
    -4.7072840136175519987205007e-02 4.4798133754572078829570003e-03
    -4.2592568383581698443496322e-02 4.4807151108790322160957409e-03
    -3.8111439955052478500974189e-02 4.4815267244567051282966297e-03
    -3.3629544980744649484982745e-02 4.4822481998660198762585161e-03
    -2.9146973606233722575709066e-02 4.4828795225957457012966323e-03
    -2.4663815990698764735178372e-02 4.4834206799479175281630816e-03
    -2.0180162305111673526347715e-02 4.4838716610380883673125929e-03
    -1.5696102730420422033397188e-02 4.4842324567955530942309217e-03
    -1.1211727455737223457798990e-02 4.4845030599635288606763517e-03
    -6.7271266765245177690624168e-03 4.4846834650993016788134149e-03
    -2.2423905927790662058474158e-03 4.4847736685743296372597122e-03
    2.2423905927790662058474158e-03 4.4847736685743296372597122e-03
    6.7271266765245177690624168e-03 4.4846834650993016788134149e-03
    1.1211727455737223457798990e-02 4.4845030599635288606763517e-03
    1.5696102730420422033397188e-02 4.4842324567955530942309217e-03
    2.0180162305111673526347715e-02 4.4838716610380883673125929e-03
    2.4663815990698764735178372e-02 4.4834206799479175281630816e-03
    2.9146973606233722575709066e-02 4.4828795225957457012966323e-03
    3.3629544980744649484982745e-02 4.4822481998660198762585161e-03
    3.8111439955052478500974189e-02 4.4815267244567051282966297e-03
    4.2592568383581698443496322e-02 4.4807151108790322160957409e-03
    4.7072840136175519987205007e-02 4.4798133754572078829570003e-03
    5.1552165099906635536974875e-02 4.4788215363280817898905539e-03
    5.6030453180892419995373643e-02 4.4777396134407830910473258e-03
    6.0507614306102325274494547e-02 4.4765676285563240494047577e-03
    6.4983558425176432549186245e-02 4.4753056052471568149186965e-03
    6.9458195512229342627463780e-02 4.4739535688966955082057630e-03
    7.3931435567667125319246679e-02 4.4725115466988200896292227e-03
    7.8403188619992875141306854e-02 4.4709795676573134415310307e-03
    8.2873364727618248348761654e-02 4.4693576625852889094847598e-03
    8.7341873980673945254693535e-02 4.4676458641045675365677248e-03
    9.1808626502814461023405102e-02 4.4658442066450127969079453e-03
    9.6273532453029789235365854e-02 4.4639527264438557882519909e-03
    1.0073650202744940940213780e-01 4.4619714615449545050407387e-03
    1.0519744546115022965082630e-01 4.4599004517980279579947300e-03
    1.0965627302996099057263990e-01 4.4577397388578582013152207e-03
    1.1411289505226755724986987e-01 4.4554893661834567980539745e-03
    1.1856722189081722595993540e-01 4.4531493790371792437787768e-03
    1.2301916395451978347708177e-01 4.4507198244838307166215685e-03
    1.2746863170025130029827665e-01 4.4482007513896954994936372e-03
    1.3191553563465446829994221e-01 4.4455922104215846168973059e-03
    1.3635978631593861920556776e-01 4.4428942540457871945847046e-03
    1.4080129435567845241905616e-01 4.4401069365270400338130408e-03
    1.4523997042061242979116287e-01 4.4372303139274199904051876e-03
    1.4967572523443883891758333e-01 4.4342644441052172718520374e-03
    1.5410846957961307768236736e-01 4.4312093867137835809244528e-03
    1.5853811429913960973081544e-01 4.4280652032003186766018210e-03
    1.6296457029836675101108767e-01 4.4248319568046439245745383e-03
    1.6738774854677854198037323e-01 4.4215097125579238060422149e-03
    1.7180756007978578714379125e-01 4.4180985372813657424684308e-03
    1.7622391600051459659148634e-01 4.4145984995848618070990454e-03
    1.8063672748159484426899724e-01 4.4110096698656278343952941e-03
    1.8504590576694612824582009e-01 4.4073321203067757426130591e-03
    1.8945136217356448038451333e-01 4.4035659248758589681682629e-03
    1.9385300809330410776176734e-01 4.3997111593233918791501225e-03
    1.9825075499466160433570394e-01 4.3957679011813292901944550e-03
    2.0264451442455458374247712e-01 4.3917362297614904662057533e-03
    2.0703419801010256029449863e-01 4.3876162261539865955262130e-03
    2.1141971746040272295275031e-01 4.3834079732255719352718160e-03
    2.1580098456830743014478458e-01 4.3791115556179828136040832e-03
    2.2017791121219731920177765e-01 4.3747270597462402028088313e-03
    2.2455040935775436028443153e-01 4.3702545737969002506706495e-03
    2.2891839105973099677271421e-01 4.3656941877262944035065217e-03
    2.3328176846372164487952716e-01 4.3610459932587122833247406e-03
    2.3764045380792672701453228e-01 4.3563100838845524725995162e-03
    2.4199435942492036888396001e-01 4.3514865548584516150021351e-03
    2.4634339774341257323087007e-01 4.3465755031973658112365300e-03
    2.5068748129000989477432881e-01 4.3415770276786129835966399e-03
    2.5502652269097653148222093e-01 4.3364912288378961585655347e-03
    2.5936043467398955941760619e-01 4.3313182089672729729867839e-03
    2.6368913006989669334245718e-01 4.3260580721130982920219488e-03
    2.6801252181446771905015680e-01 4.3207109240739416736176715e-03
    2.7233052295014642529835669e-01 4.3152768723984493218215341e-03
    2.7664304662779876098355203e-01 4.3097560263831791518840397e-03
    2.8095000610846093680450508e-01 4.3041484970704167734023571e-03
    2.8525131476508286398896530e-01 4.2984543972459246866102411e-03
    2.8954688608427203710959930e-01 4.2926738414366845397740491e-03
    2.9383663366803114413983167e-01 4.2868069459085916816931672e-03
    2.9812047123549778593343262e-01 4.2808538286641176218161320e-03
    3.0239831262468003236776326e-01 4.2748146094399265201846383e-03
    3.0667007179418792617298095e-01 4.2686894097044812690366911e-03
    3.1093566282496531982815213e-01 4.2624783526555984000672161e-03
    3.1519499992201699400951043e-01 4.2561815632179518173461297e-03
    3.1944799741613549848295861e-01 4.2497991680405956121946431e-03
    3.2369456976562277494835485e-01 4.2433312954943732536738921e-03
    3.2793463155801227948415999e-01 4.2367780756693727492456780e-03
    3.3216809751178510978775194e-01 4.2301396403722889977272459e-03
    3.3639488247808746468336949e-01 4.2234161231237809380756332e-03
    3.4061490144244083166924497e-01 4.2166076591557818606381680e-03
    3.4482806952645367326582004e-01 4.2097143854087862996360414e-03
    3.4903430198952806184919950e-01 4.2027364405290840165818622e-03
    3.5323351423056509323927799e-01 4.1956739648659861774415702e-03
    3.5742562178966352792741645e-01 4.1885271004689873450277382e-03
    3.6161054034982176297319256e-01 4.1812959910849161956902797e-03
    3.6578818573863247642918850e-01 4.1739807821550541436228521e-03
    3.6995847392997532887548573e-01 4.1665816208121863109536953e-03
    3.7412132104570755553041295e-01 4.1590986558776646409008038e-03
    3.7827664335735089462531278e-01 4.1515320378583972851793682e-03
    3.8242435728777413039836119e-01 4.1438819189438379914092181e-03
    3.8656437941287757897868005e-01 4.1361484530029147752006047e-03
    3.9069662646326763777437918e-01 4.1283317955809334387495824e-03
    3.9482101532593349979549657e-01 4.1204321038964550685812505e-03
    3.9893746304591748419454689e-01 4.1124495368381310325678335e-03
    4.0304588682798497822901140e-01 4.1043842549615102213711459e-03
    4.0714620403828893913100728e-01 4.0962364204858055238833714e-03
    4.1123833220603150939709280e-01 4.0880061972906264755600603e-03
    4.1532218902512429981044306e-01 4.0796937509127023657740096e-03
    4.1939769235584140050221436e-01 4.0712992485425281499744621e-03
    4.2346476022647383530284060e-01 4.0628228590210199028254223e-03
    4.2752331083497724373287951e-01 4.0542647528361086886605769e-03
    4.3157326255061706499205343e-01 4.0456251021193162173417157e-03
    4.3561453391561111292418218e-01 4.0369040806422940709241587e-03
    4.3964704364676626680008553e-01 4.0281018638133212969587582e-03
    4.4367071063711610579005651e-01 4.0192186286737872566443563e-03
    4.4768545395754927307407911e-01 4.0102545538946180944672726e-03
    4.5169119285843972733118790e-01 4.0012098197726945342234117e-03
    4.5568784677126900062305026e-01 3.9920846082272185006978305e-03
    4.5967533531024740156567532e-01 3.9828791027960563225773782e-03
    4.6365357827393199885435138e-01 3.9735934886320481082555567e-03
    4.6762249564683749936122581e-01 3.9642279524992815598061213e-03
    4.7158200760104723725518738e-01 3.9547826827693380313810678e-03
    4.7553203449781783485050823e-01 3.9452578694174995563304087e-03
    4.7947249688918214260979767e-01 3.9356537040189367923637143e-03
    4.8340331551954623945377421e-01 3.9259703797448440576456363e-03
    4.8732441132728337995771994e-01 3.9162080913585578870184101e-03
    4.9123570544632555456843193e-01 3.9063670352116513020956923e-03
    4.9513711920774744479345486e-01 3.8964474092399651976303954e-03
    4.9902857414135165514679215e-01 3.8864494129596405952442595e-03
    5.0290999197724384206509285e-01 3.8763732474631023249000972e-03
    5.0678129464740995224758535e-01 3.8662191154150088792662565e-03
    5.1064240428728369103339446e-01 3.8559872210481914260593594e-03
    5.1449324323731515651303425e-01 3.8456777701595255998523459e-03
    5.1833373404453042532225027e-01 3.8352909701058102327764221e-03
    5.2216379946409208212543263e-01 3.8248270297995823341352661e-03
    5.2598336246084964606950507e-01 3.8142861597049272995296754e-03
    5.2979234621089221768386324e-01 3.8036685718332353435544579e-03
    5.3359067410309157786230116e-01 3.7929744797389488251970757e-03
    5.3737826974064339946579594e-01 3.7822040985152531947233179e-03
    5.4115505694260446212240367e-01 3.7713576447897614353499396e-03
    5.4492095974542453795663732e-01 3.7604353367201568715538773e-03
    5.4867590240447439153825826e-01 3.7494373939897947776989096e-03
    5.5241980939556856178285216e-01 3.7383640378032957467258068e-03
    5.5615260541648581238405313e-01 3.7272154908820887518616782e-03
    5.5987421538848236579610784e-01 3.7159919774599385958180342e-03
    5.6358456445780258370348292e-01 3.7046937232784317266254259e-03
    5.6728357799718420739765179e-01 3.6933209555824290937220411e-03
    5.7097118160736048952941246e-01 3.6818739031155077283397059e-03
    5.7464730111855444327773057e-01 3.6703527961153554863560355e-03
    5.7831186259197264742937250e-01 3.6587578663091358671666242e-03
    5.8196479232129150194197109e-01 3.6470893469088324495563924e-03
    5.8560601683413981977111007e-01 3.6353474726065564646970696e-03
    5.8923546289357675576070505e-01 3.6235324795698235778029783e-03
    5.9285305749956462850747130e-01 3.6116446054368115878285828e-03
    5.9645872789043674622178060e-01 3.5996840893115725906747926e-03
    6.0005240154436256805325911e-01 3.5876511717592247593944244e-03
    6.0363400618080398363218819e-01 3.5755460948011198382689990e-03
    6.0720346976197092647709042e-01 3.5633691019099694372029852e-03
    6.1076072049426977095265556e-01 3.5511204380049496420745925e-03
    6.1430568682974717731326564e-01 3.5388003494467757011066045e-03
    6.1783829746752882972060661e-01 3.5264090840327480882998579e-03
    6.2135848135525451052529888e-01 3.5139468909917625713545952e-03
    6.2486616769050606912117019e-01 3.5014140209793051007614029e-03
    6.2836128592223272626426933e-01 3.4888107260724063665713324e-03
    6.3184376575216827376380024e-01 3.4761372597645729363990963e-03
    6.3531353713624749701693872e-01 3.4633938769606828315950064e-03
    6.3877053028601249451412514e-01 3.4505808339718719961186633e-03
    6.4221467567001810916593740e-01 3.4376983885103609174527239e-03
    6.4564590401522947704648914e-01 3.4247467996842899211340505e-03
    6.4906414630841602342314900e-01 3.4117263279924928826014607e-03
    6.5246933379753857540350737e-01 3.3986372353192679032773960e-03
    6.5586139799313358800247897e-01 3.3854797849291002817539997e-03
    6.5924027066968915455902334e-01 3.3722542414613772450426854e-03
    6.6260588386701835261760607e-01 3.3589608709250575770133285e-03
    6.6595816989162581744920999e-01 3.3455999406933256343221750e-03
    6.6929706131806887547952556e-01 3.3321717194982150046789204e-03
    6.7262249099031456989195021e-01 3.3186764774251957359207932e-03
    6.7593439202308958080323009e-01 3.3051144859077511567457819e-03
    6.7923269780322581556930572e-01 3.2914860177219108956780946e-03
    6.8251734199100011490912721e-01 3.2777913469807713232884172e-03
    6.8578825852146896302485857e-01 3.2640307491289782641785688e-03
    6.8904538160579675842853931e-01 3.2502045009371819533905423e-03
    6.9228864573257897774283265e-01 3.2363128804964850539216403e-03
    6.9551798566916078758737285e-01 3.2223561672128338620457733e-03
    6.9873333646294777388163766e-01 3.2083346418014073442304124e-03
    7.0193463344271289638953704e-01 3.1942485862809654080518573e-03
    7.0512181221989789214887878e-01 3.1800982839681828616418269e-03
    7.0829480868990712938426668e-01 3.1658840194719460765793162e-03
    7.1145355903339813075092479e-01 3.1516060786876305688242184e-03
    7.1459799971756421399504688e-01 3.1372647487913521251179372e-03
    7.1772806749741302478895477e-01 3.1228603182341862706805635e-03
    7.2084369941703840822810889e-01 3.1083930767363695223115538e-03
    7.2394483281088617410148345e-01 3.0938633152814754880000070e-03
    7.2703140530501531024754058e-01 3.0792713261105541036610767e-03
    7.3010335481835231252745189e-01 3.0646174027162604615315278e-03
    7.3316061956393907550477707e-01 3.0499018398369498114575116e-03
    7.3620313805017678632225397e-01 3.0351249334507495770962837e-03
    7.3923084908206193599511380e-01 3.0202869807696049175849140e-03
    7.4224369176241744572308789e-01 3.0053882802333021714846772e-03
    7.4524160549311724288656933e-01 2.9904291315034658461924355e-03
    7.4822452997630561899455870e-01 2.9754098354575278517764048e-03
    7.5119240521560926016064741e-01 2.9603306941826876275136726e-03
    7.5414517151734494770920492e-01 2.9451920109698193593616988e-03
    7.5708276949171915415348622e-01 2.9299940903073887384089335e-03
    7.6000514005402330930394328e-01 2.9147372378753129071315975e-03
    7.6291222442582184992687644e-01 2.8994217605388225740548069e-03
    7.6580396413613460726566018e-01 2.8840479663422885665824413e-03
    7.6868030102261264424612364e-01 2.8686161645030232303366091e-03
    7.7154117723270820849990059e-01 2.8531266654050636139006158e-03
    7.7438653522483846813884156e-01 2.8375797805929251632628407e-03
    7.7721631776954280823588306e-01 2.8219758227653406711110673e-03
    7.8003046795063346596776910e-01 2.8063151057689645316572946e-03
    7.8282892916634128077646437e-01 2.7905979445920674544834306e-03
    7.8561164513045289581327779e-01 2.7748246553581956165557276e-03
    7.8837855987344362951318999e-01 2.7589955553198146354088038e-03
    7.9112961774360357480873063e-01 2.7431109628519253530731792e-03
    7.9386476340815570473807838e-01 2.7271711974456679106193491e-03
    7.9658394185437009227257477e-01 2.7111765797018880924662554e-03
    7.9928709839066958142694830e-01 2.6951274313246884918593427e-03
    8.0197417864772979623211313e-01 2.6790240751149692680077585e-03
    8.0464512857957359859284452e-01 2.6628668349639194636024886e-03
    8.0729989446465699742816469e-01 2.6466560358465143938666575e-03
    8.0993842290695006180811788e-01 2.6303920038149761727319653e-03
    8.1256066083701139479700259e-01 2.6140750659922112539290673e-03
    8.1516655551305472471312896e-01 2.5977055505652392984605736e-03
    8.1775605452201027834036040e-01 2.5812837867785843118384470e-03
    8.2032910578057882666769274e-01 2.5648101049276579750657579e-03
    8.2288565753627918031298805e-01 2.5482848363521134853193484e-03
    8.2542565836848913463086319e-01 2.5317083134291798809933915e-03
    8.2794905718947975348243062e-01 2.5150808695669811879125355e-03
    8.3045580324544276962228651e-01 2.4984028391978295446929081e-03
    8.3294584611751176783656092e-01 2.4816745577714927409318513e-03
    8.3541913572277581856440065e-01 2.4648963617484556838654886e-03
    8.3787562231528678324821158e-01 2.4480685885931484715993811e-03
    8.4031525648706084652417303e-01 2.4311915767671579864661613e-03
    8.4273798916907072253934530e-01 2.4142656657224238758718826e-03
    8.4514377163223375344358601e-01 2.3972911958944072112476231e-03
    8.4753255548839157018647938e-01 2.3802685086952457028941499e-03
    8.4990429269128375810993248e-01 2.3631979465068889655066098e-03
    8.5225893553751330689038923e-01 2.3460798526742051607618667e-03
    8.5459643666750717549973615e-01 2.3289145714980807010119346e-03
    8.5691674906646797538201099e-01 2.3117024482284974015722234e-03
    8.5921982606531999149268586e-01 2.2944438290575784079872168e-03
    8.6150562134164809791059270e-01 2.2771390611126388937857090e-03
    8.6377408892062879086637395e-01 2.2597884924491899207021905e-03
    8.6602518317595578167811254e-01 2.2423924720439488041112686e-03
    8.6825885883075637483585751e-01 2.2249513497878173860777817e-03
    8.7047507095850373826095847e-01 2.2074654764788399254060725e-03
    8.7267377498391918155817848e-01 2.1899352038151527477527480e-03
    8.7485492668386966030880103e-01 2.1723608843879096096107784e-03
    8.7701848218825650960184248e-01 2.1547428716741862456118817e-03
    8.7916439798089818236093151e-01 2.1370815200298736401263167e-03
    8.8129263090040554917692361e-01 2.1193771846825552526705216e-03
    8.8340313814105009271315794e-01 2.1016302217243491151643653e-03
    8.8549587725362433054954181e-01 2.0838409881047633732953361e-03
    8.8757080614629635384460471e-01 2.0660098416235002198593218e-03
    8.8962788308545570625796017e-01 2.0481371409232698027613750e-03
    8.9166706669655348971303965e-01 2.0302232454825746427173883e-03
    8.9368831596493392144253676e-01 2.0122685156084758363592702e-03
    8.9569159023665978480721606e-01 1.9942733124293514530844806e-03
    8.9767684921932955344203720e-01 1.9762379978876267426490809e-03
    8.9964405298288829815334111e-01 1.9581629347325008733138318e-03
    9.0159316196043082225486387e-01 1.9400484865126482996594559e-03
    9.0352413694899691432027566e-01 1.9218950175689025663267051e-03
    9.0543693911036049470908438e-01 1.9037028930269346907455663e-03
    9.0733152997181054644215692e-01 1.8854724787899042239697200e-03
    9.0920787142692471860527803e-01 1.8672041415311048905001368e-03
    9.1106592573633582432535150e-01 1.8488982486865846409374026e-03
    9.1290565552849134434154621e-01 1.8305551684477561636549270e-03
    9.1472702380040482594836249e-01 1.8121752697539956870886879e-03
    9.1652999391839951037752598e-01 1.7937589222852166285365749e-03
    9.1831452961884640906475852e-01 1.7753064964544380387873046e-03
    9.2008059500889272097623461e-01 1.7568183634003315795463207e-03
    9.2182815456718381064149526e-01 1.7382948949797596104038799e-03
    9.2355717314457808075900402e-01 1.7197364637602950612066399e-03
    9.2526761596485385119592593e-01 1.7011434430127274266414394e-03
    9.2695944862540868847133879e-01 1.6825162067035522809460568e-03
    9.2863263709795118572287720e-01 1.6638551294874545042473679e-03
    9.3028714772918563724601881e-01 1.6451605866997726437817029e-03
    9.3192294724148894147219835e-01 1.6264329543489433257952292e-03
    9.3354000273357984340805160e-01 1.6076726091089474021678107e-03
    9.3513828168118051653578959e-01 1.5888799283117253056951679e-03
    9.3671775193767115030851755e-01 1.5700552899395952242966867e-03
    9.3827838173473598892826431e-01 1.5511990726176452543710882e-03
    9.3982013968300293083046881e-01 1.5323116556061207836625382e-03
    9.4134299477267413536196727e-01 1.5133934187927951774133017e-03
    9.4284691637415030118773984e-01 1.4944447426853291888138031e-03
    9.4433187423864661802497267e-01 1.4754660084036210453140026e-03
    9.4579783849880105783825002e-01 1.4564575976721339835295854e-03
    9.4724477966927489447357402e-01 1.4374198928122235673077167e-03
    9.4867266864734633990963175e-01 1.4183532767344489596006429e-03
    9.5008147671349552076947020e-01 1.3992581329308674976258375e-03
    9.5147117553198179429330139e-01 1.3801348454673253648983255e-03
    9.5284173715141451399546213e-01 1.3609837989757285307834689e-03
    9.5419313400531458047026945e-01 1.3418053786463130395284482e-03
    9.5552533891266877574821592e-01 1.3225999702198892784416051e-03
    9.5683832507847732529171481e-01 1.3033679599800925344443847e-03
    9.5813206609429157900592600e-01 1.2841097347456062287285317e-03
    9.5940653593874647420136625e-01 1.2648256818623868861367621e-03
    9.6066170897808289552699534e-01 1.2455161891958719741491102e-03
    9.6189755996666392867666673e-01 1.2261816451231773335284192e-03
    9.6311406404748201026677634e-01 1.2068224385252889711145352e-03
    9.6431119675265963842036854e-01 1.1874389587792388064263482e-03
    9.6548893400394086850013764e-01 1.1680315957502765151360125e-03
    9.6664725211317581443637437e-01 1.1486007397840248926701445e-03
    9.6778612778279737849373987e-01 1.1291467816986319652045045e-03
    9.6890553810628932129844770e-01 1.1096701127769103070730417e-03
    9.7000546056864733746039064e-01 1.0901711247584676845595597e-03
    9.7108587304683191554488531e-01 1.0706502098318268576676582e-03
    9.7214675381021331546094189e-01 1.0511077606265416956032865e-03
    9.7318808152100833019915171e-01 1.0315441702052981134268839e-03
    9.7420983523470983111991472e-01 1.0119598320560106718968285e-03
    9.7521199440050809759128470e-01 9.9235514008391115421681139e-04
    9.7619453886170381995412981e-01 9.7273048860362402393692216e-04
    9.7715744885611366399302824e-01 9.5308627233123993129082496e-04
    9.7810070501646784180138638e-01 9.3342288637637761856968854e-04
    9.7902428837079980006308233e-01 9.1374072623424116342677470e-04
    9.7992818034282747063912211e-01 8.9404018777766638016030187e-04
    9.8081236275232719368233347e-01 8.7432166724916331786848778e-04
    9.8167681781549931407937493e-01 8.5458556125295192813090539e-04
    9.8252152814532589530926998e-01 8.3483226674698667382767958e-04
    9.8334647675192077276307145e-01 8.1506218103498189690020448e-04
    9.8415164704287072527932878e-01 7.9527570175842797039339471e-04
    9.8493702282356965227450019e-01 7.7547322688860127054555349e-04
    9.8570258829754442420068017e-01 7.5565515471857100469194046e-04
    9.8644832806677229530833984e-01 7.3582188385520343708312563e-04
    9.8717422713199121098170963e-01 7.1597381321115885054878758e-04
    9.8788027089300134431226752e-01 6.9611134199689035130520498e-04
    9.8856644514895852804414744e-01 6.7623486971263996325731682e-04
    9.8923273609866091415909750e-01 6.5634479614043493914121719e-04
    9.8987913034082575247651903e-01 6.3644152133608743269338470e-04
    9.9050561487436006302687019e-01 6.1652544562120089287338143e-04
    9.9111217709862164948475538e-01 5.9659696957518036121437266e-04
    9.9169880481367378433077420e-01 5.7665649402726456163731372e-04
    9.9226548622053001302845132e-01 5.5670442004856959122838234e-04
    9.9281220992139318504143830e-01 5.3674114894416698022527878e-04
    9.9333896491988438182119125e-01 5.1676708224520316852651503e-04
    9.9384574062126618265722300e-01 4.9678262170107492703280405e-04
    9.9433252683265593852013353e-01 4.7678816927169146095213947e-04
    9.9479931376323316172261002e-01 4.5678412711985843163917109e-04
    9.9524609202443881095234701e-01 4.3677089760383425815978842e-04
    9.9567285263016636065458442e-01 4.1674888327014565159514548e-04
    9.9607958699694720827721994e-01 3.9671848684678294536670728e-04
    9.9646628694412919813316876e-01 3.7668011123697140796542926e-04
    9.9683294469404826187997060e-01 3.5663415951382569256053467e-04
    9.9717955287219695037492784e-01 3.3658103491638908383309925e-04
    9.9750610450738885770505249e-01 3.1652114084790224394702629e-04
    9.9781259303192360032852548e-01 2.9645488087773699554389184e-04
    9.9809901228175612608595202e-01 2.7638265874956383354871248e-04
    9.9836535649667978997712225e-01 2.5630487840044917065632224e-04
    9.9861162032053651937957284e-01 2.3622194399989265996796239e-04
    9.9883779880148382268600926e-01 2.1613426002684216816612794e-04
    9.9904388739237171002116611e-01 1.9604223142283443213733374e-04
    9.9922988195134376798733911e-01 1.7594626390712142906605497e-04
    9.9939577874291152248531489e-01 1.5584676466203183911249375e-04
    9.9954157444010394151234777e-01 1.3574414394147891441307074e-04
    9.9966726612928136219693442e-01 1.1563881924509692518312359e-04
    9.9977285132240067966336028e-01 9.5531227691402493119691985e-05
    9.9985832799377194479717446e-01 7.5421869998298820669119236e-05
    9.9992369471822484250367324e-01 5.5311513823237843051703472e-05
    9.9996895141230646153474027e-01 3.5202627719042746520707737e-05
    9.9999410721789827594108147e-01 1.5122766974076512616403586e-05];

X=xw(:,1); W=xw(:,2);









% function save_moments(moms)
%
% %--------------------------------------------------------------------------
% % OBJECT:
% %--------------------------------------------------------------------------
% % Save moments on a file.
% %--------------------------------------------------------------------------
% % INPUT:
% %--------------------------------------------------------------------------
% % moms: column vector of moments.
% %--------------------------------------------------------------------------
%
% CC=clock;
% clock_str=strcat(num2str(CC(4)),'_',...
%     num2str(CC(5)),'_',...
%     num2str(floor(CC(6))));
% filestr=strcat('moms_saved_GG',clock_str,'.m');
% fid = fopen(filestr,'w');
% fprintf(fid,'function moms=\n');
% fprintf(fid,filestr);
% fprintf(fid,'moms=[ \n');
% fprintf(fid,'%1.25e \n',moms);
% fprintf(fid,']; \n');
% fclose(fid);