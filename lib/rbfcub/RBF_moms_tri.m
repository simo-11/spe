
function moms=RBF_moms_tri(P,centers,RBF_type,RBF_scale)

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
% 13: phi=@(r) exp(-r);                     % Matern beta_1=3/2.
% 14: phi=@(r) (1+r).*exp(-r);              % Matern beta_2=5/2.
% 15: phi=@(r) (max(0,(1-r)))^10*(429*r^4 + 450*r^3 + 210*r^2 + 50*r + 5)
%                                           % Wendland W8
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
%% Date: November 2019-October 10, 2020.
%--------------------------------------------------------------------------

if nargin < 1, error('Domain not defined.'); end
if nargin < 2, error('Centers not defined.'); end
if nargin < 3, warning('RBF not declared. Using TPS.'); RBF_type=5; end
if nargin < 4
    warning('RBF_scale not declared. Using RBF_scale=1.');
    RBF_scale=1;
end

if length(RBF_scale) < size(centers,1)
    RBF_scale=RBF_scale*ones(size(centers,1),1);
end

tri = triangulation(P);
TCL=tri.ConnectivityList;
X = tri.Points(:,1);
Y = tri.Points(:,2);

M=size(centers,1);
N=size(TCL,1);

for jj=1:M
    
    C=centers(jj,:);
    RBF_scaleL=RBF_scale(jj);
    for ii=1:N
        TTloc=TCL(ii,:);
        Xloc=X(TTloc');
        Yloc=Y(TTloc');
        momsL(ii)=RBFcub_ABC([Xloc Yloc],C,RBF_type,RBF_scaleL);
    end
    moms(jj,1)=sum(momsL);
    
end

moms=real(moms);





%--------------------------------------------------------------------------
% RBFcub_ABC
%--------------------------------------------------------------------------

function I=RBFcub_ABC(vertices,center,RBF_type,RBF_scale)

%--------------------------------------------------------------------------
% OBJECT.
% Integral over the general triangle having vertices "vertices",
% using as integrand the Wendland RBF with center "center" and scale
% "RBF_scale".
%--------------------------------------------------------------------------
% INPUT
% vertices: triangle vertices as 3 x 2 matrix (no need to repeat the first
% vertex!);
% center: RBF center;
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
%
% RBF_scale: a positive scalar: all the RBFs have the same scale.
%            a vector of the same number of rows of centers: scale of the
%            RBF for each center.
%--------------------------------------------------------------------------
% OUTPUT
% I: Integral over the general triangle having vertices "vertices",
% using as integrand the Wendland RBF with center "center" and scale
% "RBF_scale"
%--------------------------------------------------------------------------
% RELATED WORK.
%--------------------------------------------------------------------------
% Alvise Sommariva, Marco Vianello
% "Numerical cubature on scattered datas by Radial Basis Functions",
% Computing 76 (2005), 295-310.
%--------------------------------------------------------------------------
% Copyrights.
%--------------------------------------------------------------------------
%% Copyright (C) 2006-2020 Alvise Sommariva, Marco Vianello.
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
%% Date: May 2004 - October 14, 2020.
%--------------------------------------------------------------------------

% SETTING DEFAULTS.

if nargin < 3, RBF_scale=1; end

if polyarea(vertices(:,1),vertices(:,2)) < 10^(-20), I=0; return; end

happy=inpolygon(center(1),center(2),vertices(:,1),vertices(:,2));

if happy == 1 % the triangle contains the center
    %fprintf('\n \t ..... CENTER IN TRIANGLE .... \n');
    I_OAB=RBFcub_OAB(center,vertices(1,:),vertices(2,:),RBF_type,RBF_scale);
    I_OBC=RBFcub_OAB(center,vertices(2,:),vertices(3,:),RBF_type,RBF_scale);
    I_OCA=RBFcub_OAB(center,vertices(3,:),vertices(1,:),RBF_type,RBF_scale);
    I=I_OAB+I_OBC+I_OCA;
else % the triangle does not contain the center
    [see,no_see,P]=see_all_sides(center,vertices);
    
    if length(no_see) == 0 % see all vertices
        
        imin=[];
        for k=1:2
            % could do with setdiff, but it is faster
            if k==1
                iother=[2 3];
            else
                iother=[1 3];
            end
            
            X=vertices(k,1); Y=vertices(k,2);
            
            XV=[center(1); vertices(iother,1)];
            YV=[center(2); vertices(iother,2)];
            
            in=inpolygon(X,Y,XV,YV);
            if in == 1
                imin=k; break;
            end
            
        end
        
        if length(imin) == 0
            imin=3; iother=[1 2];
        end
        
        
        A=vertices(imin,:);
        B=vertices(iother(1),:);
        C=vertices(iother(2),:);
        
        I_OBC=RBFcub_OAB(center,B,C,RBF_type,RBF_scale);
        I_OCA=RBFcub_OAB(center,A,C,RBF_type,RBF_scale);
        I_OAB=RBFcub_OAB(center,A,B,RBF_type,RBF_scale);
        
        I=I_OBC-I_OCA-I_OAB;
        
    else
        % ... one obscured vertex, say A ...
        A=vertices(no_see,:);
        B=vertices(see(1),:);
        C=vertices(see(2),:);
        
        I_OCA=RBFcub_OAB(center,A,C,RBF_type,RBF_scale);
        I_OAB=RBFcub_OAB(center,A,B,RBF_type,RBF_scale);
        I_OBC=RBFcub_OAB(center,B,C,RBF_type,RBF_scale);
        
        I=I_OCA+I_OAB-I_OBC;
    end
    
end




%--------------------------------------------------------------------------
% "see_all_sides" and subroutines.
%--------------------------------------------------------------------------
function [see,no_see,P]=see_all_sides(C,vertices)

V1=vertices(1,:);
V2=vertices(2,:);
V3=vertices(3,:);

see=[1 2 3]; no_see=[];

[happyL,P]=test_vertex_see(C,V1,V2,V3);
if (happyL == 0)
    see=[2 3]; no_see=1;
    return;
end

[happyL,P]=test_vertex_see(C,V2,V1,V3);
if (happyL == 0)
    see=[1 3]; no_see=2;
    return;
end

[happyL,P]=test_vertex_see(C,V3,V1,V2);
if (happyL == 0)
    see=[1 2]; no_see=3;
    return;
end

P=[];


function [happy,P]=test_vertex_see(V1,V2,V3,V4)

A=[(V2-V1)' -(V4-V3)'];
b=(V3-V1)';
th=A\b;
P=V3+th(2)*(V4-V3);
happy= ~((th(1) >= 0 & th(1) <= 1) & (th(2) >= 0 & th(2) <= 1));








%--------------------------------------------------------------------------
% "RBFcub_OAB" and subroutines.
%--------------------------------------------------------------------------

function I=RBFcub_OAB(O,A,B,RBF_type,RBF_scale)

%--------------------------------------------------------------------------
% OBJECT.
% Integral over the right triangle OAB, using Wendland RBF with
% scale "RBF_scale".
%--------------------------------------------------------------------------
% INPUT
% 0: vertex of the triangle OAB and RBF center.
% H: vertex of the triangle OAB, with angle OHB equal to "pi/2".
% B vertex of the triangle OAB.
% RBF_type:
% RBF_scale:
%--------------------------------------------------------------------------
% OUTPUT
% I: integral of the Wendland RBF with center in "O" that may not be the
% origin (0,0) and scale "RBF_scale" over the right triangle OAB.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% RELATED WORK.
%--------------------------------------------------------------------------
% Alvise Sommariva, Marco Vianello
% "Numerical cubature on scattered datas by Radial Basis Functions",
% Computing 76 (2005), 295-310.
%--------------------------------------------------------------------------
% Copyrights.
%--------------------------------------------------------------------------
%% Copyright (C) 2006-2019 Alvise Sommariva, Marco Vianello.
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
%% Date: May 2004 - October 25, 2019.
%--------------------------------------------------------------------------

% SETTING DEFAULTS.

if nargin < 4, RBF_scale=1; end

misAO=norm(O-A);
if misAO == 0, I=0; return; end
misBO=norm(O-B);
if misBO == 0, I=0; return; end

[H,d] = find_projection_to_line(O, A, B);
happy = projection_in_segment(H, A, B);

misHO=norm(O-H);
misHA=norm(A-H);
misHB=norm(B-H);

I_OHA=RBFcub_OHB(misHO,misAO,misHA,RBF_type,RBF_scale);
I_OHB=RBFcub_OHB(misHO,misBO,misHB,RBF_type,RBF_scale);

if happy == 1
    I=I_OHA+I_OHB;
else
    misHA=norm(H-A);
    misHB=norm(H-B);
    if misHA <= misHB
        I=I_OHB-I_OHA;
    else
        I=I_OHA-I_OHB;
    end
end




function res = projection_in_segment(pt, V1, V2)

sx=pt-V2;
dx=V1-V2;

if abs(dx(1)) > 0
    theta=sx(1)/dx(1);
else
    if abs(dx(2)) > 0
        theta=sx(2)/dx(2);
    else % V1=V2
        if P == V2
            theta=0;
        else
            theta=+inf;
        end
    end
end

res=(theta >= 0) & (theta <= 1);



function [P,d] = find_projection_to_line(C, V1, V2)

% % Compute distance
% a = v1 - v2;
% b = pt - v2;
% d = norm(cross(a,b)) / norm(a);
% % Normalize couple of vectors
% ae = a / norm(a);
% be = b / norm(b);
% % Two cross products give direction of perpendicular
% h = cross(ae,cross(ae,be));
% he = h / norm(h);
% % Perpendicular is its base vector times length
% perp = he*d;

if size(C,1) < size(C,2), C=C'; flag=0; end
if size(V1,1) < size(V1,2), V1=V1'; end
if size(V2,1) < size(V2,2), V2=V2'; end

u=C-V1;
v=(V2-V1)/norm(V2-V1);
s=u'*v;
P=V1+s*v;
d=norm(P-C);

if flag == 0, P=P'; end





















function I=RBFcub_OHB(misHO,misBO,misHB,RBF_type,RBF_scale)

%--------------------------------------------------------------------------
% Input:
%--------------------------------------------------------------------------
% RBF_type integer corresponding to RBF
%    1: phi=@(r) (1+r.*r).^(1/2);             % Multiquadric
%    2: phi=@(r) exp(-r.*r);                  % Gaussian
%    3: phi=@(r) (1+r.*r).^(-1/2);            % Inverse Multiquadric
%    4: phi=@(r) (1+4*r).*(max(0,(1-r)).^4);  % Wendland 2
%    5: phi=@(r) r.*r.*log(r+eps*(1-sign(r))); % TPS
%    6: phi=@(r) r.^3;                        % polyharmonic spline
%    7: phi=@(r) r.^5;                        % polyharmonic spline
%    8: phi=@(r) r.^7;                        % polyharmonic spline
%    9: phi=@(r) max(0,(1-r)).^2;             % Wendland W0
%    10: phi=@(r) (35*r.^2+18*r+3).*(max(0,(1-r)))^6;         % Wendland W4
%    11: phi=@(r) (32*r.^3+25*r.^2+8*r+1).*(max(0,(1-r).)^8;  % Wendland W6
%    13 phi=@(r) exp(-r);                     % Matern beta_1=3/2.
%    14 phi=@(r) (1+r).*exp(-r);              % Matern beta_2=5/2.
%    15: phi=@(r) (max(0,(1-r)))^10*(429*r^4 + 450*r^3 + 210*r^2 + 50*r+5)
%                                           % Wendland W8
%--------------------------------------------------------------------------
% Output:
%--------------------------------------------------------------------------
% I: integral on a right triangle OHA.
%--------------------------------------------------------------------------

tol=10^(-16);

if (misHO <= tol) | (misBO <= tol) | (misHB <= tol)
    I=0; return;
end


switch RBF_type
    case 1
        I=MQ_integral(misBO,misHB,misHO,RBF_scale);     % Multiquadric
    case 2
        I=gaussian_integral(misBO,misHB,misHO,RBF_scale); % Gaussian
    case 3
        I=IMQ_integral(misBO,misHB,misHO,RBF_scale); % Inverse Multiquadric
    case 4
        I=W2_integral_19(misBO,misHB,misHO,RBF_scale);  % Wendland 2
    case 5
        I=TPS_integral(misBO,misHB,misHO,RBF_scale); % TPS
    case 6
        I=R3_integral(misBO,misHB,misHO,RBF_scale);   % polyharmonic spline
    case 7
        I=R5_integral(misBO,misHB,misHO,RBF_scale);   % polyharmonic spline
    case 8
        I=R7_integral(misBO,misHB,misHO,RBF_scale);   % polyharmonic spline
    case 9
        I=W_integral_19(misBO,misHB,misHO,RBF_scale,0);     % Wendland W0
    case 10
        I=W_integral_19(misBO,misHB,misHO,RBF_scale,4);  % Wendland W4
    case 11
        I=W_integral_19(misBO,misHB,misHO,RBF_scale,6);  % Wendland W6
    case 12
        I=MW0_integral_19(misBO,misHB,misHO,RBF_scale); % Missing Wendland
    case 13
        I=M1_integral_19(misBO,misHB,misHO,RBF_scale); % Matern beta_1=3/2.
    case 14
        I=M2_integral_19(misBO,misHB,misHO,RBF_scale); % Missing Wendland
    case 15
        I=W_integral_19(misBO,misHB,misHO,RBF_scale,8);  % Wendland W8
    otherwise
        error('RBF type not implemented')
end











function I=MQ_integral(misOA,misHA,misHO,RBF_scale)

%--------------------------------------------------------------------------
% Object.
%--------------------------------------------------------------------------
% It computes the integral of the Multiquadric RBF with center the origin,
% RBF scale equal to "RBF_scale", on the triangle OHA, right in H.
%--------------------------------------------------------------------------

method=1;

switch method
    
    case 1
        
        % ....................... some assignments ........................
        
        
        
        a=misOA; b=misHA; c=misHO; sig=RBF_scale;
        bs=b/sig; cs=c/sig; gam=c/b;
        
        % .................................................................
        % Note:
        % The primitive is defined as the sum of two functions, i.e. "fI"
        % and "fII", that must be evaluated in "bs" and "0".
        % .................................................................
        
        % .................... first function analysis ....................
        
        % The code below computes, just faster:
        %
        % fI=@(c,t) (1/2)*(...
        %     (2/3)*c*t*sqrt(c^2+t^2+1)+(1/3)*(t^2+3)*t*log(sqrt(c^2+t^2+1)+c) + ...
        %     (-1/6)*c*(c^2-3)*log(sqrt(c^2+t^2+1)+t)+...
        %     (1/2)*c*(c^2+1)*log(sqrt(c^2+t^2+1)+t) + ...
        %     (-i/3)*log( (6*(-i*c^2+t-i))/(c^2*(t-i)) - ...
        %     (6*i*sqrt(c^2+t^2+1))/(c*(t-i))  ) + ...
        %     (i/3)*log( (6*i*sqrt(c^2+t^2+1))/(c*(t+i)) + ...
        %     (6*(i*c^2+t+i))/(c^2*(t+i))   ) + ...
        %     -t^3/9-2*t/3+(2/3)*atan(t)+0*c);
        %
        % term(1)=(fI(cs,bs)-fI(cs,0));
        
        c=cs; t=bs;
        termL(1)=(1/2)*(...
            (2/3)*c*t*sqrt(c^2+t^2+1)+ ...
            (1/3)*(t^2+3)*t*log(sqrt(c^2+t^2+1)+c) + ...
            (-1/6)*c*(c^2-3)*log(sqrt(c^2+t^2+1)+t)+...
            (1/2)*c*(c^2+1)*log(sqrt(c^2+t^2+1)+t) + ...
            (-i/3)*log( (6*(-i*c^2+t-i))/(c^2*(t-i)) - ...
            (6*i*sqrt(c^2+t^2+1))/(c*(t-i))  ) + ...
            (i/3)*log( (6*i*sqrt(c^2+t^2+1))/(c*(t+i)) + ...
            (6*(i*c^2+t+i))/(c^2*(t+i))   ) + ...
            -t^3/9-2*t/3+(2/3)*atan(t)+0*c);
        
        c=cs; t=0;
        termL(2)=(1/2)*( (-1/6)*c*(c^2-3)*log(sqrt(c^2+1))+...
            (1/2)*c*(c^2+1)*log(sqrt(c^2+1)) + ...
            (-i/3)*log( (6*(-i*c^2-i))/(c^2*(-i)) - ...
            (6*i*sqrt(c^2+1))/(c*(-i))  ) + ...
            (i/3)*log( (6*i*sqrt(c^2+1))/(c*(i)) + ...
            (6*(i*c^2+i))/(c^2*(i))   ) );
        
        
        
        term(1)=termL(1)-termL(2);
        
        
        % .................... second function analysis ...................
        
        
        % The code below computes, just faster:
        %
        % fII=@(gam,t) (1/18)*(...
        %     3*gam*t^2*sqrt( (gam^2+1)*t^2+1 )+...
        %     3*(t^2+3)*t*log(sqrt( (gam^2+1)*t^2+1 )+gam*t)+...
        %     3*i*log(6*(-gam*sqrt( (gam^2+1)*t^2+1 )+gam^2*t+t-i) /...
        %     (gam^2*(t-i)) )+...
        %     3*i*log(-6*(gam*sqrt( (gam^2+1)*t^2+1 )+gam^2*t+t+i) /...
        %     (gam^2*(t+i)) )+...
        %     -t^3-6*t+6*atan(t) );
        % term(2)=(fII(gam,bs)-fII(gam,0));
        
        c=gam; t=bs;
        termL(1)=(1/18)*(...
            3*gam*t^2*sqrt( (gam^2+1)*t^2+1 )+...
            3*(t^2+3)*t*log(sqrt((gam^2+1)*t^2+1 )+gam*t)+...
            3*i*log(6*(-gam*sqrt((gam^2+1)*t^2+1 )+gam^2*t+t-i)/(gam^2*(t-i)) )+...
            3*i*log(-6*(gam*sqrt((gam^2+1)*t^2+1 )+gam^2*t+t+i)/(gam^2*(t+i)) )+...
            -t^3-6*t+6*atan(t) );
        
        c=gam; t=0;
        termL(2)=(1/18)*(3*i*log(6*(-gam-i)/(gam^2*(-i)))+3*i*...
            log(-6*(gam+i)/(gam^2*(i))));
        
        term(2)=termL(1)-termL(2);
        
        
        % .................... integral on triangle OHA ...........................
        
        I=sig^2*real(term(1)-term(2));
        
        
    case 2
        
        r0=misHO/RBF_scale; r1=misOA/RBF_scale; theta=acos(r0/r1);
        
        f=@(t,c) (sec(t)*(c^2 + (cos(t))^2)*sqrt(c^2*(sec(t))^2 + 1)*...
            (sqrt(2)*c^2*(c^2 + 3)*(cos(t))* ...
            atanh((sqrt(2)*sqrt(c^2)*sin(t))/sqrt(2*c^2 + cos(2*t) + 1)) + ...
            sqrt(c^2)*(c^2*sin(t)*sqrt(2*c^2 + cos(2*t) + 1) - ...
            2*i*sqrt(2)*(cos(t))^2* log(sqrt(2*c^2 + cos(2*t) + 1) +...
            i*sqrt(2)*sin(t)))))/(3*sqrt(c^2)*(2*c^2 + cos(2*t) + 1)^(3/2));
        
        I=RBF_scale^2*(f(theta,r0)-f(0,r0));
        
end











function I=gaussian_integral(misOA,misHA,misHO,RBF_scale)

%--------------------------------------------------------------------------
% Object.
%--------------------------------------------------------------------------
% It computes the integral of the gaussian RBF with center the origin, RBF
% scale equal to "RBF_scale", on the triangle OHA, right in H.
%--------------------------------------------------------------------------

a=misOA; b=misHA; c=misHO; sig=RBF_scale;

IL(1)=(pi/4)*(sig^2)*erf(c/sig)*erf(b/sig);

value=primitive2eval(b,c,sig);
IL(2)=(-sqrt(pi)*sig^2/2)*value;

I=sum(IL);


function value=primitive2eval(b,c,sig)

C=c/b; B=b/sig;

% Note:
% Formula provided by using Mathematica command:
%            integral (exp(-s^2)*erf(c*s)) ds
% at Wolfram alpha integrator site.

% TB=tfn(sqrt(2)*abs(C)*B, 1/abs(C) );
%IB=(sqrt(pi)/(2*abs(C)))*(4*C*TB+abs(C)*erf(B)*erf(C*B)+C);
TB=TFNG(sqrt(2)*C*B, 1/C );
IB=2*sqrt(pi)*TB+(sqrt(pi)/2)*erf(B)*erf(C*B)+(sqrt(pi)/2);

T0=TFNG(0, 1/C );
I0=2*sqrt(pi)*T0+(sqrt(pi)/2);

value=IB-I0;


function value = TFNG ( x,a )

XW=[3.138447481065978e-03     6.276874353199929e-03
    9.415218790093920e-03     6.276627047804665e-03
    1.569161914468320e-02     6.276132446757833e-03
    2.196740125811324e-02     6.275390569546442e-03
    2.824231786802290e-02     6.274401445400047e-03
    3.451612174615035e-02     6.273165113289592e-03
    4.078856570807577e-02     6.271681621925886e-03
    4.705940262295932e-02     6.269951029757671e-03
    5.332838542327661e-02     6.267973404969322e-03
    5.959526711455566e-02     6.265748825478173e-03
    6.585980078510630e-02     6.263277378931428e-03
    7.212173961574986e-02     6.260559162702727e-03
    7.838083688954177e-02     6.257594283888281e-03
    8.463684600149458e-02     6.254382859302696e-03
    9.088952046829173e-02     6.250925015474328e-03
    9.713861393800019e-02     6.247220888640323e-03
    1.033838801997763e-01     6.243270624741241e-03
    1.096250731935659e-01     6.239074379415299e-03
    1.158619470198002e-01     6.234632317992259e-03
    1.220942559490822e-01     6.229944615486894e-03
    1.283217544318700e-01     6.225011456592109e-03
    1.345441971081505e-01     6.219833035671644e-03
    1.407613388171065e-01     6.214409556752433e-03
    1.469729346067760e-01     6.208741233516565e-03
    1.531787397437033e-01     6.202828289292851e-03
    1.593785097225811e-01     6.196670957048045e-03
    1.655720002758839e-01     6.190269479377645e-03
    1.717589673834924e-01     6.183624108496350e-03
    1.779391672823075e-01     6.176735106228118e-03
    1.841123564758542e-01     6.169602743995839e-03
    1.902782917438753e-01     6.162227302810670e-03
    1.964367301519146e-01     6.154609073260926e-03
    2.025874290608877e-01     6.146748355500671e-03
    2.087301461366425e-01     6.138645459237846e-03
    2.148646393595063e-01     6.130300703722111e-03
    2.209906670338228e-01     6.121714417732241e-03
    2.271079877974721e-01     6.112886939563174e-03
    2.332163606313832e-01     6.103818617012701e-03
    2.393155448690276e-01     6.094509807367736e-03
    2.454053002059024e-01     6.084960877390250e-03
    2.514853867089998e-01     6.075172203302837e-03
    2.575555648262564e-01     6.065144170773865e-03
    2.636155953959952e-01     6.054877174902302e-03
    2.696652396563472e-01     6.044371620202133e-03
    2.757042592546568e-01     6.033627920586431e-03
    2.817324162568760e-01     6.022646499351047e-03
    2.877494731569358e-01     6.011427789157935e-03
    2.937551928861057e-01     5.999972232018107e-03
    2.997493388223336e-01     5.988280279274203e-03
    3.057316747995684e-01     5.976352391582735e-03
    3.117019651170646e-01     5.964189038895907e-03
    3.176599745486687e-01     5.951790700443126e-03
    3.236054683520887e-01     5.939157864712104e-03
    3.295382122781403e-01     5.926291029429615e-03
    3.354579725799772e-01     5.913190701541888e-03
    3.413645160223022e-01     5.899857397194632e-03
    3.472576098905530e-01     5.886291641712705e-03
    3.531370220000754e-01     5.872493969579401e-03
    3.590025207052670e-01     5.858464924415415e-03
    3.648538749087072e-01     5.844205058957404e-03
    3.706908540702600e-01     5.829714935036222e-03
    3.765132282161589e-01     5.814995123554770e-03
    3.823207679480669e-01     5.800046204465528e-03
    3.881132444521138e-01     5.784868766747680e-03
    3.938904295079139e-01     5.769463408383918e-03
    3.996520954975543e-01     5.753830736336884e-03
    4.053980154145667e-01     5.737971366525251e-03
    4.111279628728678e-01     5.721885923799469e-03
    4.168417121156816e-01     5.705575041917127e-03
    4.225390380244313e-01     5.689039363517990e-03
    4.282197161276125e-01     5.672279540098694e-03
    4.338835226096336e-01     5.655296231987053e-03
    4.395302343196360e-01     5.638090108316064e-03
    4.451596287802860e-01     5.620661846997532e-03
    4.507714841965398e-01     5.603012134695359e-03
    4.563655794643815e-01     5.585141666798498e-03
    4.619416941795375e-01     5.567051147393552e-03
    4.674996086461555e-01     5.548741289237027e-03
    4.730391038854647e-01     5.530212813727262e-03
    4.785599616444008e-01     5.511466450875990e-03
    4.840619644042067e-01     5.492502939279596e-03
    4.895448953890008e-01     5.473323026090000e-03
    4.950085385743206e-01     5.453927466985228e-03
    5.004526786956304e-01     5.434317026139632e-03
    5.058771012568061e-01     5.414492476193794e-03
    5.112815925385832e-01     5.394454598224067e-03
    5.166659396069802e-01     5.374204181711820e-03
    5.220299303216853e-01     5.353742024512319e-03
    5.273733533444160e-01     5.333068932823296e-03
    5.326959981472462e-01     5.312185721153188e-03
    5.379976550208994e-01     5.291093212289042e-03
    5.432781150830127e-01     5.269792237264101e-03
    5.485371702863653e-01     5.248283635325061e-03
    5.537746134270766e-01     5.226568253898995e-03
    5.589902381527688e-01     5.204646948559990e-03
    5.641838389706982e-01     5.182520582995407e-03
    5.693552112558506e-01     5.160190028971876e-03
    5.745041512590039e-01     5.137656166300933e-03
    5.796304561147554e-01     5.114919882804373e-03
    5.847339238495148e-01     5.091982074279251e-03
    5.898143533894621e-01     5.068843644462598e-03
    5.948715445684692e-01     5.045505504995817e-03
    5.999052981359867e-01     5.021968575388761e-03
    6.049154157648942e-01     4.998233782983507e-03
    6.099017000593141e-01     4.974302062917807e-03
    6.148639545623893e-01     4.950174358088268e-03
    6.198019837640234e-01     4.925851619113181e-03
    6.247155931085830e-01     4.901334804295079e-03
    6.296045890025642e-01     4.876624879582970e-03
    6.344687788222193e-01     4.851722818534293e-03
    6.393079709211461e-01     4.826629602276550e-03
    6.441219746378392e-01     4.801346219468651e-03
    6.489106003032011e-01     4.775873666261967e-03
    6.536736592480160e-01     4.750212946261077e-03
    6.584109638103827e-01     4.724365070484226e-03
    6.631223273431084e-01     4.698331057323498e-03
    6.678075642210626e-01     4.672111932504684e-03
    6.724664898484908e-01     4.645708729046874e-03
    6.770989206662871e-01     4.619122487221756e-03
    6.817046741592265e-01     4.592354254512627e-03
    6.862835688631563e-01     4.565405085573125e-03
    6.908354243721447e-01     4.538276042185677e-03
    6.953600613455899e-01     4.510968193219663e-03
    6.998573015152848e-01     4.483482614589307e-03
    7.043269676924420e-01     4.455820389211280e-03
    7.087688837746735e-01     4.427982606962043e-03
    7.131828747529304e-01     4.399970364634893e-03
    7.175687667183971e-01     4.371784765896767e-03
    7.219263868693439e-01     4.343426921244745e-03
    7.262555635179350e-01     4.314897947962303e-03
    7.305561260969925e-01     4.286198970075286e-03
    7.348279051667178e-01     4.257331118307632e-03
    7.390707324213666e-01     4.228295530036816e-03
    7.432844406958796e-01     4.199093349249037e-03
    7.474688639724699e-01     4.169725726494153e-03
    7.516238373871631e-01     4.140193818840337e-03
    7.557491972362932e-01     4.110498789828503e-03
    7.598447809829525e-01     4.080641809426451e-03
    7.639104272633952e-01     4.050624053982780e-03
    7.679459758933953e-01     4.020446706180541e-03
    7.719512678745577e-01     3.990110954990627e-03
    7.759261454005825e-01     3.959617995624947e-03
    7.798704518634827e-01     3.928969029489322e-03
    7.837840318597541e-01     3.898165264136148e-03
    7.876667311964989e-01     3.867207913216835e-03
    7.915183968974998e-01     3.836098196433975e-03
    7.953388772092481e-01     3.804837339493290e-03
    7.991280216069219e-01     3.773426574055344e-03
    8.028856808003169e-01     3.741867137687014e-03
    8.066117067397294e-01     3.710160273812727e-03
    8.103059526217874e-01     3.678307231665474e-03
    8.139682728952365e-01     3.646309266237590e-03
    8.175985232666736e-01     3.614167638231310e-03
    8.211965607062316e-01     3.581883614009088e-03
    8.247622434532156e-01     3.549458465543722e-03
    8.282954310216878e-01     3.516893470368221e-03
    8.317959842060022e-01     3.484189911525477e-03
    8.352637650862903e-01     3.451349077517717e-03
    8.386986370338935e-01     3.418372262255733e-03
    8.421004647167474e-01     3.385260765007906e-03
    8.454691141047138e-01     3.352015890349011e-03
    8.488044524748606e-01     3.318638948108821e-03
    8.521063484166915e-01     3.285131253320500e-03
    8.553746718373237e-01     3.251494126168792e-03
    8.586092939666133e-01     3.217728891938000e-03
    8.618100873622285e-01     3.183836880959783e-03
    8.649769259146711e-01     3.149819428560733e-03
    8.681096848522449e-01     3.115677875009760e-03
    8.712082407459723e-01     3.081413565465301e-03
    8.742724715144563e-01     3.047027849922302e-03
    8.773022564286910e-01     3.012522083159052e-03
    8.802974761168184e-01     2.977897624683778e-03
    8.832580125688314e-01     2.943155838681110e-03
    8.861837491412231e-01     2.908298093958310e-03
    8.890745705615831e-01     2.873325763891355e-03
    8.919303629331384e-01     2.838240226370824e-03
    8.947510137392416e-01     2.803042863747606e-03
    8.975364118478033e-01     2.767735062778438e-03
    9.002864475156713e-01     2.732318214571274e-03
    9.030010123929542e-01     2.696793714530465e-03
    9.056799995272898e-01     2.661162962301793e-03
    9.083233033680599e-01     2.625427361717316e-03
    9.109308197705484e-01     2.589588320740065e-03
    9.135024460000440e-01     2.553647251408565e-03
    9.160380807358894e-01     2.517605569781210e-03
    9.185376240754719e-01     2.481464695880463e-03
    9.210009775381599e-01     2.445226053636912e-03
    9.234280440691834e-01     2.408891070833172e-03
    9.258187280434574e-01     2.372461179047622e-03
    9.281729352693497e-01     2.335937813598016e-03
    9.304905729923920e-01     2.299322413484917e-03
    9.327715498989343e-01     2.262616421335017e-03
    9.350157761197427e-01     2.225821283344289e-03
    9.372231632335403e-01     2.188938449221011e-03
    9.393936242704907e-01     2.151969372128651e-03
    9.415270737156246e-01     2.114915508628618e-03
    9.436234275122095e-01     2.077778318622864e-03
    9.456826030650609e-01     2.040559265296379e-03
    9.477045192437970e-01     2.003259815059539e-03
    9.496890963860347e-01     1.965881437490332e-03
    9.516362563005288e-01     1.928425605276460e-03
    9.535459222702526e-01     1.890893794157321e-03
    9.554180190554201e-01     1.853287482865870e-03
    9.572524728964510e-01     1.815608153070356e-03
    9.590492115168765e-01     1.777857289315959e-03
    9.608081641261872e-01     1.740036378966301e-03
    9.625292614226219e-01     1.702146912144848e-03
    9.642124355958984e-01     1.664190381676218e-03
    9.658576203298855e-01     1.626168283027366e-03
    9.674647508052152e-01     1.588082114248681e-03
    9.690337637018371e-01     1.549933375914974e-03
    9.705645972015134e-01     1.511723571066379e-03
    9.720571909902542e-01     1.473454205149153e-03
    9.735114862606945e-01     1.435126785956388e-03
    9.749274257144106e-01     1.396742823568639e-03
    9.763049535641789e-01     1.358303830294469e-03
    9.776440155361730e-01     1.319811320610912e-03
    9.789445588721034e-01     1.281266811103870e-03
    9.802065323312958e-01     1.242671820408439e-03
    9.814298861927109e-01     1.204027869149173e-03
    9.826145722569035e-01     1.165336479880304e-03
    9.837605438479229e-01     1.126599177025921e-03
    9.848677558151522e-01     1.087817486820118e-03
    9.859361645350886e-01     1.048992937247145e-03
    9.869657279130641e-01     1.010127057981583e-03
    9.879564053849051e-01     9.712213803285645e-04
    9.889081579185339e-01     9.322774371641146e-04
    9.898209480155092e-01     8.932967628756748e-04
    9.906947397125084e-01     8.542808933029134e-04
    9.915294985827499e-01     8.152313656789921e-04
    9.923251917373580e-01     7.761497185725336e-04
    9.930817878266696e-01     7.370374918306677e-04
    9.937992570414846e-01     6.978962265237658e-04
    9.944775711142633e-01     6.587274648928495e-04
    9.951167033202715e-01     6.195327503013223e-04
    9.957166284786823e-01     5.803136271938493e-04
    9.962773229536414e-01     5.410716410674078e-04
    9.967987646553128e-01     5.018083384637189e-04
    9.972809330409360e-01     4.625252670007075e-04
    9.977238091159455e-01     4.232239754783270e-04
    9.981273754352638e-01     3.839060141334781e-04
    9.984916161049857e-01     3.445729352122161e-04
    9.988165167849478e-01     3.052262942671818e-04
    9.991020646933555e-01     2.658676532633462e-04
    9.993482486165864e-01     2.264985887095048e-04
    9.995550589335447e-01     1.871207158506811e-04
    9.997224876879449e-01     1.477357747755901e-04
    9.998505288592006e-01     1.083460286909200e-04
    9.999391798145371e-01     6.895707282668985e-05
    9.999884567522129e-01     2.962364448548283e-05];
X=XW(:,1); W=XW(:,2);
fX=exp(-0.5*x^2*(1+a^2*X.^2))./(1+a^2*X.^2);
value=(a/(2*pi))*(W'*fX);






function value = TFNC ( x,a )

f=@(t) (1/(2*pi))*(exp(-x^2*(1+t.^2)/2))./(1+t.^2);
F=chebfun(f,[0 a],'splitting','on');
value=sum(F);



function value = tfn ( x, fx )

%**************************************************************************
%
%% TFN calculates the T-function of Owen.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    20 January 2008
%
%  Author:
%
%    Original FORTRAN77 version by JC Young, Christoph Minder.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    MA Porter, DJ Winstanley,
%    Remark AS R30:
%    A Remark on Algorithm AS76:
%    An Integral Useful in Calculating Noncentral T and Bivariate
%    Normal Probabilities,
%    Applied Statistics,
%    Volume 28, Number 1, 1979, page 113.
%
%    JC Young, Christoph Minder,
%    Algorithm AS 76:
%    An Algorithm Useful in Calculating Non-Central T and
%    Bivariate Normal Distributions,
%    Applied Statistics,
%    Volume 23, Number 3, 1974, pages 455-457.
%
%  Parameters:
%
%    Input, real X, FX, the parameters of the function.
%
%    Output, real VALUE, the value of the T-function.
%
ng = 5;

r = [ ...
    0.1477621, ...
    0.1346334, ...
    0.1095432, ...
    0.0747257, ...
    0.0333357 ];
tp = 0.159155;
tv1 = 1.0E-35;
tv2 = 15.0;
tv3 = 15.0;
tv4 = 1.0E-05;
u = [ ...
    0.0744372, ...
    0.2166977, ...
    0.3397048, ...
    0.4325317, ...
    0.4869533 ];
%
%  Test for X near zero.
%
if ( abs ( x ) < tv1 )
    value = tp * atan ( fx );
    return
end
%
%  Test for large values of abs(X).
%
if ( tv2 < abs ( x ) )
    value = 0.0;
    return
end
%
%  Test for FX near zero.
%
if ( abs ( fx ) < tv1 )
    value = 0.0;
    return
end
%
%  Test whether abs ( FX ) is so large that it must be truncated.
%
xs = - 0.5 * x * x;
x2 = fx;
fxs = fx * fx;
%
%  Computation of truncation point by Newton iteration.
%
if ( tv3 <= log ( 1.0 + fxs ) - xs * fxs )
    
    x1 = 0.5 * fx;
    fxs = 0.25 * fxs;
    
    while ( 1 )
        
        rt = fxs + 1.0;
        
        x2 = x1 + ( xs * fxs + tv3 - log ( rt ) ) ...
            / ( 2.0 * x1 * ( 1.0 / rt - xs ) );
        
        fxs = x2 * x2;
        
        if ( abs ( x2 - x1 ) < tv4 )
            break
        end
        
        x1 = x2;
        
    end
    
end
%
%  Gaussian quadrature.
%
a=fx;
XW=[3.138447481065978e-03     6.276874353199929e-03
    9.415218790093920e-03     6.276627047804665e-03
    1.569161914468320e-02     6.276132446757833e-03
    2.196740125811324e-02     6.275390569546442e-03
    2.824231786802290e-02     6.274401445400047e-03
    3.451612174615035e-02     6.273165113289592e-03
    4.078856570807577e-02     6.271681621925886e-03
    4.705940262295932e-02     6.269951029757671e-03
    5.332838542327661e-02     6.267973404969322e-03
    5.959526711455566e-02     6.265748825478173e-03
    6.585980078510630e-02     6.263277378931428e-03
    7.212173961574986e-02     6.260559162702727e-03
    7.838083688954177e-02     6.257594283888281e-03
    8.463684600149458e-02     6.254382859302696e-03
    9.088952046829173e-02     6.250925015474328e-03
    9.713861393800019e-02     6.247220888640323e-03
    1.033838801997763e-01     6.243270624741241e-03
    1.096250731935659e-01     6.239074379415299e-03
    1.158619470198002e-01     6.234632317992259e-03
    1.220942559490822e-01     6.229944615486894e-03
    1.283217544318700e-01     6.225011456592109e-03
    1.345441971081505e-01     6.219833035671644e-03
    1.407613388171065e-01     6.214409556752433e-03
    1.469729346067760e-01     6.208741233516565e-03
    1.531787397437033e-01     6.202828289292851e-03
    1.593785097225811e-01     6.196670957048045e-03
    1.655720002758839e-01     6.190269479377645e-03
    1.717589673834924e-01     6.183624108496350e-03
    1.779391672823075e-01     6.176735106228118e-03
    1.841123564758542e-01     6.169602743995839e-03
    1.902782917438753e-01     6.162227302810670e-03
    1.964367301519146e-01     6.154609073260926e-03
    2.025874290608877e-01     6.146748355500671e-03
    2.087301461366425e-01     6.138645459237846e-03
    2.148646393595063e-01     6.130300703722111e-03
    2.209906670338228e-01     6.121714417732241e-03
    2.271079877974721e-01     6.112886939563174e-03
    2.332163606313832e-01     6.103818617012701e-03
    2.393155448690276e-01     6.094509807367736e-03
    2.454053002059024e-01     6.084960877390250e-03
    2.514853867089998e-01     6.075172203302837e-03
    2.575555648262564e-01     6.065144170773865e-03
    2.636155953959952e-01     6.054877174902302e-03
    2.696652396563472e-01     6.044371620202133e-03
    2.757042592546568e-01     6.033627920586431e-03
    2.817324162568760e-01     6.022646499351047e-03
    2.877494731569358e-01     6.011427789157935e-03
    2.937551928861057e-01     5.999972232018107e-03
    2.997493388223336e-01     5.988280279274203e-03
    3.057316747995684e-01     5.976352391582735e-03
    3.117019651170646e-01     5.964189038895907e-03
    3.176599745486687e-01     5.951790700443126e-03
    3.236054683520887e-01     5.939157864712104e-03
    3.295382122781403e-01     5.926291029429615e-03
    3.354579725799772e-01     5.913190701541888e-03
    3.413645160223022e-01     5.899857397194632e-03
    3.472576098905530e-01     5.886291641712705e-03
    3.531370220000754e-01     5.872493969579401e-03
    3.590025207052670e-01     5.858464924415415e-03
    3.648538749087072e-01     5.844205058957404e-03
    3.706908540702600e-01     5.829714935036222e-03
    3.765132282161589e-01     5.814995123554770e-03
    3.823207679480669e-01     5.800046204465528e-03
    3.881132444521138e-01     5.784868766747680e-03
    3.938904295079139e-01     5.769463408383918e-03
    3.996520954975543e-01     5.753830736336884e-03
    4.053980154145667e-01     5.737971366525251e-03
    4.111279628728678e-01     5.721885923799469e-03
    4.168417121156816e-01     5.705575041917127e-03
    4.225390380244313e-01     5.689039363517990e-03
    4.282197161276125e-01     5.672279540098694e-03
    4.338835226096336e-01     5.655296231987053e-03
    4.395302343196360e-01     5.638090108316064e-03
    4.451596287802860e-01     5.620661846997532e-03
    4.507714841965398e-01     5.603012134695359e-03
    4.563655794643815e-01     5.585141666798498e-03
    4.619416941795375e-01     5.567051147393552e-03
    4.674996086461555e-01     5.548741289237027e-03
    4.730391038854647e-01     5.530212813727262e-03
    4.785599616444008e-01     5.511466450875990e-03
    4.840619644042067e-01     5.492502939279596e-03
    4.895448953890008e-01     5.473323026090000e-03
    4.950085385743206e-01     5.453927466985228e-03
    5.004526786956304e-01     5.434317026139632e-03
    5.058771012568061e-01     5.414492476193794e-03
    5.112815925385832e-01     5.394454598224067e-03
    5.166659396069802e-01     5.374204181711820e-03
    5.220299303216853e-01     5.353742024512319e-03
    5.273733533444160e-01     5.333068932823296e-03
    5.326959981472462e-01     5.312185721153188e-03
    5.379976550208994e-01     5.291093212289042e-03
    5.432781150830127e-01     5.269792237264101e-03
    5.485371702863653e-01     5.248283635325061e-03
    5.537746134270766e-01     5.226568253898995e-03
    5.589902381527688e-01     5.204646948559990e-03
    5.641838389706982e-01     5.182520582995407e-03
    5.693552112558506e-01     5.160190028971876e-03
    5.745041512590039e-01     5.137656166300933e-03
    5.796304561147554e-01     5.114919882804373e-03
    5.847339238495148e-01     5.091982074279251e-03
    5.898143533894621e-01     5.068843644462598e-03
    5.948715445684692e-01     5.045505504995817e-03
    5.999052981359867e-01     5.021968575388761e-03
    6.049154157648942e-01     4.998233782983507e-03
    6.099017000593141e-01     4.974302062917807e-03
    6.148639545623893e-01     4.950174358088268e-03
    6.198019837640234e-01     4.925851619113181e-03
    6.247155931085830e-01     4.901334804295079e-03
    6.296045890025642e-01     4.876624879582970e-03
    6.344687788222193e-01     4.851722818534293e-03
    6.393079709211461e-01     4.826629602276550e-03
    6.441219746378392e-01     4.801346219468651e-03
    6.489106003032011e-01     4.775873666261967e-03
    6.536736592480160e-01     4.750212946261077e-03
    6.584109638103827e-01     4.724365070484226e-03
    6.631223273431084e-01     4.698331057323498e-03
    6.678075642210626e-01     4.672111932504684e-03
    6.724664898484908e-01     4.645708729046874e-03
    6.770989206662871e-01     4.619122487221756e-03
    6.817046741592265e-01     4.592354254512627e-03
    6.862835688631563e-01     4.565405085573125e-03
    6.908354243721447e-01     4.538276042185677e-03
    6.953600613455899e-01     4.510968193219663e-03
    6.998573015152848e-01     4.483482614589307e-03
    7.043269676924420e-01     4.455820389211280e-03
    7.087688837746735e-01     4.427982606962043e-03
    7.131828747529304e-01     4.399970364634893e-03
    7.175687667183971e-01     4.371784765896767e-03
    7.219263868693439e-01     4.343426921244745e-03
    7.262555635179350e-01     4.314897947962303e-03
    7.305561260969925e-01     4.286198970075286e-03
    7.348279051667178e-01     4.257331118307632e-03
    7.390707324213666e-01     4.228295530036816e-03
    7.432844406958796e-01     4.199093349249037e-03
    7.474688639724699e-01     4.169725726494153e-03
    7.516238373871631e-01     4.140193818840337e-03
    7.557491972362932e-01     4.110498789828503e-03
    7.598447809829525e-01     4.080641809426451e-03
    7.639104272633952e-01     4.050624053982780e-03
    7.679459758933953e-01     4.020446706180541e-03
    7.719512678745577e-01     3.990110954990627e-03
    7.759261454005825e-01     3.959617995624947e-03
    7.798704518634827e-01     3.928969029489322e-03
    7.837840318597541e-01     3.898165264136148e-03
    7.876667311964989e-01     3.867207913216835e-03
    7.915183968974998e-01     3.836098196433975e-03
    7.953388772092481e-01     3.804837339493290e-03
    7.991280216069219e-01     3.773426574055344e-03
    8.028856808003169e-01     3.741867137687014e-03
    8.066117067397294e-01     3.710160273812727e-03
    8.103059526217874e-01     3.678307231665474e-03
    8.139682728952365e-01     3.646309266237590e-03
    8.175985232666736e-01     3.614167638231310e-03
    8.211965607062316e-01     3.581883614009088e-03
    8.247622434532156e-01     3.549458465543722e-03
    8.282954310216878e-01     3.516893470368221e-03
    8.317959842060022e-01     3.484189911525477e-03
    8.352637650862903e-01     3.451349077517717e-03
    8.386986370338935e-01     3.418372262255733e-03
    8.421004647167474e-01     3.385260765007906e-03
    8.454691141047138e-01     3.352015890349011e-03
    8.488044524748606e-01     3.318638948108821e-03
    8.521063484166915e-01     3.285131253320500e-03
    8.553746718373237e-01     3.251494126168792e-03
    8.586092939666133e-01     3.217728891938000e-03
    8.618100873622285e-01     3.183836880959783e-03
    8.649769259146711e-01     3.149819428560733e-03
    8.681096848522449e-01     3.115677875009760e-03
    8.712082407459723e-01     3.081413565465301e-03
    8.742724715144563e-01     3.047027849922302e-03
    8.773022564286910e-01     3.012522083159052e-03
    8.802974761168184e-01     2.977897624683778e-03
    8.832580125688314e-01     2.943155838681110e-03
    8.861837491412231e-01     2.908298093958310e-03
    8.890745705615831e-01     2.873325763891355e-03
    8.919303629331384e-01     2.838240226370824e-03
    8.947510137392416e-01     2.803042863747606e-03
    8.975364118478033e-01     2.767735062778438e-03
    9.002864475156713e-01     2.732318214571274e-03
    9.030010123929542e-01     2.696793714530465e-03
    9.056799995272898e-01     2.661162962301793e-03
    9.083233033680599e-01     2.625427361717316e-03
    9.109308197705484e-01     2.589588320740065e-03
    9.135024460000440e-01     2.553647251408565e-03
    9.160380807358894e-01     2.517605569781210e-03
    9.185376240754719e-01     2.481464695880463e-03
    9.210009775381599e-01     2.445226053636912e-03
    9.234280440691834e-01     2.408891070833172e-03
    9.258187280434574e-01     2.372461179047622e-03
    9.281729352693497e-01     2.335937813598016e-03
    9.304905729923920e-01     2.299322413484917e-03
    9.327715498989343e-01     2.262616421335017e-03
    9.350157761197427e-01     2.225821283344289e-03
    9.372231632335403e-01     2.188938449221011e-03
    9.393936242704907e-01     2.151969372128651e-03
    9.415270737156246e-01     2.114915508628618e-03
    9.436234275122095e-01     2.077778318622864e-03
    9.456826030650609e-01     2.040559265296379e-03
    9.477045192437970e-01     2.003259815059539e-03
    9.496890963860347e-01     1.965881437490332e-03
    9.516362563005288e-01     1.928425605276460e-03
    9.535459222702526e-01     1.890893794157321e-03
    9.554180190554201e-01     1.853287482865870e-03
    9.572524728964510e-01     1.815608153070356e-03
    9.590492115168765e-01     1.777857289315959e-03
    9.608081641261872e-01     1.740036378966301e-03
    9.625292614226219e-01     1.702146912144848e-03
    9.642124355958984e-01     1.664190381676218e-03
    9.658576203298855e-01     1.626168283027366e-03
    9.674647508052152e-01     1.588082114248681e-03
    9.690337637018371e-01     1.549933375914974e-03
    9.705645972015134e-01     1.511723571066379e-03
    9.720571909902542e-01     1.473454205149153e-03
    9.735114862606945e-01     1.435126785956388e-03
    9.749274257144106e-01     1.396742823568639e-03
    9.763049535641789e-01     1.358303830294469e-03
    9.776440155361730e-01     1.319811320610912e-03
    9.789445588721034e-01     1.281266811103870e-03
    9.802065323312958e-01     1.242671820408439e-03
    9.814298861927109e-01     1.204027869149173e-03
    9.826145722569035e-01     1.165336479880304e-03
    9.837605438479229e-01     1.126599177025921e-03
    9.848677558151522e-01     1.087817486820118e-03
    9.859361645350886e-01     1.048992937247145e-03
    9.869657279130641e-01     1.010127057981583e-03
    9.879564053849051e-01     9.712213803285645e-04
    9.889081579185339e-01     9.322774371641146e-04
    9.898209480155092e-01     8.932967628756748e-04
    9.906947397125084e-01     8.542808933029134e-04
    9.915294985827499e-01     8.152313656789921e-04
    9.923251917373580e-01     7.761497185725336e-04
    9.930817878266696e-01     7.370374918306677e-04
    9.937992570414846e-01     6.978962265237658e-04
    9.944775711142633e-01     6.587274648928495e-04
    9.951167033202715e-01     6.195327503013223e-04
    9.957166284786823e-01     5.803136271938493e-04
    9.962773229536414e-01     5.410716410674078e-04
    9.967987646553128e-01     5.018083384637189e-04
    9.972809330409360e-01     4.625252670007075e-04
    9.977238091159455e-01     4.232239754783270e-04
    9.981273754352638e-01     3.839060141334781e-04
    9.984916161049857e-01     3.445729352122161e-04
    9.988165167849478e-01     3.052262942671818e-04
    9.991020646933555e-01     2.658676532633462e-04
    9.993482486165864e-01     2.264985887095048e-04
    9.995550589335447e-01     1.871207158506811e-04
    9.997224876879449e-01     1.477357747755901e-04
    9.998505288592006e-01     1.083460286909200e-04
    9.999391798145371e-01     6.895707282668985e-05
    9.999884567522129e-01     2.962364448548283e-05];
X=XW(:,1); W=XW(:,2);
fX=exp(-0.5*x^2*(1+a^2*X.^2))./(1+a^2*X.^2);
value=(a/(2*pi))*(W'*fX);











function I=IMQ_integral(misOA,misHA,misHO,RBF_scale)

%--------------------------------------------------------------------------
% Object.
%--------------------------------------------------------------------------
% It computes the integral of the IMQ RBF, i.e. phi(r)=(1+r*r)^(-1/2),
% with center in O, RBF scale equal to "RBF_scale", on the triangle OHA,
% right in H.
%--------------------------------------------------------------------------

method=1;

switch method
    case 1 % method from cartesian coordinates
        
        % ....................... some assignments ........................
        
        a=misOA; b=misHA; c=misHO; sig=RBF_scale;
        bs=b/sig; cs=c/sig; gam=c/b;
        
        % .................................................................
        % Note:
        % The primitive is defined as the sum of two functions, i.e. "fI"
        % and "fII", that must be evaluated in "bs" and "0".
        % .................................................................
        
        % .................... first function analysis ....................
        
        % the code below performs (quickly)
        % fI=@(c,t) t*log( sqrt(c^2+t^2+1) + c)+...
        %    c*log(sqrt(c^2+t^2+1) + t)+...
        %    (-i/2)*log((-i*c*sqrt(c^2+t^2+1)-i*c^2+t-i )/(c^2*(t-i)))+...
        %    (i/2)*log((i*c*sqrt(c^2+t^2+1)+i*c^2+t+i)/(c^2*(t+i)))+ ...
        %    (-1)*t+atan(t);
        %
        % term(1)=(fI(cs,bs)-fI(cs,0));
        
        t=bs; c=cs;
        termL(1)=t*log( sqrt(c^2+t^2+1) + c)+...
            c*log(sqrt(c^2+t^2+1) + t)+...
            (-i/2)*log((-i*c*sqrt(c^2+t^2+1)-i*c^2+t-i )/(c^2*(t-i)))+...
            (i/2)*log((i*c*sqrt(c^2+t^2+1)+i*c^2+t+i)/(c^2*(t+i)))+ ...
            (-1)*t+atan(t);
        
        t=0; c=cs;
        %         termL(2)=t*log( sqrt(c^2+t^2+1) + c)+...
        %             c*log(sqrt(c^2+t^2+1) + t)+...
        %             (-i/2)*log((-i*c*sqrt(c^2+t^2+1)-i*c^2+t-i )/(c^2*(t-i)))+...
        %             (i/2)*log((i*c*sqrt(c^2+t^2+1)+i*c^2+t+i)/(c^2*(t+i)))+ ...
        %             (-1)*t+atan(t);
        
        termL(2)=c*log(sqrt(c^2+1))+...
            (-i/2)*log((-i*c*sqrt(c^2+1)-i*c^2-i )/(c^2*(-i)))+...
            (i/2)*log((i*c*sqrt(c^2+1)+i*c^2+i)/(c^2*(i)));
        
        term(1)=termL(1)-termL(2);
        
        % .................... second function analysis ...................
        
        % the code below performs (quickly)
        %
        % fII=@(gam,t) t*log( sqrt( (gam^2+1)*t^2+1 ) + gam*t ) + ...
        %    (i/2)*log((4*(-gam*sqrt((gam^2+1)*t^2+1)+...
        %    gam^2*t+t-i))/(gam^2*(t-i)))+...
        %    (i/2)*log(-(4*(gam*sqrt((gam^2+1)*t^2+1)+...
        %    gam^2*t+t+i))/(gam^2*(t+i)))+(-1)*t+atan(t);
        % term(2)=(fII(gam,bs)-fII(gam,0));
        
        t=bs; c=gam;
        termL(1)= t*log( sqrt( (gam^2+1)*t^2+1 ) + gam*t ) + ...
            (i/2)*log((4*(-gam*sqrt((gam^2+1)*t^2+1)+...
            gam^2*t+t-i))/(gam^2*(t-i)))+...
            (i/2)*log(-(4*(gam*sqrt((gam^2+1)*t^2+1)+...
            gam^2*t+t+i))/(gam^2*(t+i)))+(-1)*t+atan(t);
        
        t=0; c=gam;
        termL(2)=(i/2)*log((4*(-gam*sqrt((gam^2+1)*t^2+1)+...
            gam^2*t+t-i))/(gam^2*(t-i)))+...
            (i/2)*log(-(4*(gam*sqrt((gam^2+1)*t^2+1)+...
            gam^2*t+t+i))/(gam^2*(t+i)));
        
        term(2)=termL(1)-termL(2);
        % .................... integral on triangle OHA ...................
        
        I=sig^2*real(term(1)-term(2));
        
    case 2 % method from polar coordinates (not working)
        
        r0=misHO/RBF_scale; r1=misOA/RBF_scale; theta=acos(r0/r1);
        
        % the code below performs (more quickly!):
        %
        % f=@(t,c) (cos(t)*sqrt(2*c^2*(sec(t))^2 + 2)* ...
        %     (sqrt(c^2 + 1)* asin(sin(t)/sqrt(c^2 + 1))* sqrt(2* c^2 + ...
        %     cos(2*t) + 1)* sqrt((2* c^2 + cos(2* t) + 1)/(c^2 + 1)) + ...
        %     sqrt(c^2)* (2* c^2 + cos(2* t) + 1)* atanh((sqrt(2)* sqrt(c^2)* ...
        %     sin(t))/sqrt(2 *c^2 + cos(2* t) + 1))))/(2* c^2 + cos(2* t) + 1)^(3/2);
        % I=f(theta,r0)-f(0,r0);
        % Note that "f(0,r0)=0".
        
        c=r0; t=theta;
        
        I=(cos(t)*sqrt(2*c^2*(sec(t))^2 + 2)* ...
            (sqrt(c^2 + 1)* asin(sin(t)/sqrt(c^2 + 1))* sqrt(2* c^2 + ...
            cos(2*t) + 1)* sqrt((2* c^2 + cos(2* t) + 1)/(c^2 + 1)) + ...
            sqrt(c^2)* (2* c^2 + cos(2* t) + 1)* atanh((sqrt(2)* sqrt(c^2)* ...
            sin(t))/sqrt(2 *c^2 + cos(2* t) + 1))))/(2* c^2 + cos(2* t) + 1)^(3/2);
        
        I=RBF_scale^2*I;
        
end











function I=W2_integral_19(misOA,misHA,misHO,RBF_scale)

%--------------------------------------------------------------------------
% Object.
%--------------------------------------------------------------------------
% It computes the integral of the Wendland RBF (W2) with center the origin,
% RBF scale equal to "RBF_scale", on the triangle OHA, right in H.
%--------------------------------------------------------------------------

% ....................... some assignments ................................

a=misOA; b=misHA; c=misHO; sig=RBF_scale;
bs=b/sig; cs=c/sig; gam=c/b;

mu0=pi/7; % measure in R^2 of the Wendland RBF with scale equal to 1.

% ..........  the region is union of a sector and a right triangle ........

if RBF_scale < a & RBF_scale > c
    
    % ... computing measure of sector ...
    
    thetaM=asin(c/sig);
    thetam=asin(c/a);
    mu_sector=sig^2*mu0*(thetaM-thetam)/(2*pi);
    
    % ... computing RBF integral on a triangle, in the support ...
    
    % we reduce the computation to a right triangle, OHJ, right in H,
    % with "mis0J=sig".
    misOJ=RBF_scale;
    misHJ=sqrt(RBF_scale^2-c^2);
    int_insupp=W2_integral_19_insupp(misOJ,misHJ,misHO,RBF_scale);
    %int_insupp =W_integral_19_insupp(misOJ,misHJ,misHO,RBF_scale,2);
    
    I=mu_sector+int_insupp;
    return;
    
end

% ..........  the region is a sector ..........

% ... in this instance the radius of the sector is RBF_scale ...
if c >= RBF_scale
    
    thetaM=pi/2;
    thetam=asin(c/a);
    I=sig^2*mu0*(thetaM-thetam)/(2*pi);
    
    
else
    % ..........  the region is a right triangle ..........
    I=W2_integral_19_insupp(misOA,misHA,misHO,RBF_scale);
    % I=W_integral_19_insupp(misOA,misHA,misHO,RBF_scale,2);
end





function I=W2_integral_19_insupp(misOA,misHA,misHO,RBF_scale)

%--------------------------------------------------------------------------
% Object.
%--------------------------------------------------------------------------
% It computes the integral of the Wendland W2 RBF with center the origin O,
% where the RBF scale equal to "RBF_scale", on the triangle OHA, right in H.
% Here it is supposed that the triangle is inside the support.
%--------------------------------------------------------------------------

r0=misHO/RBF_scale; r1=misOA/RBF_scale; theta=acos(r0/r1);

% the code below performs (more quickly!):
%
%         W2_int=@(a,t) 2/21*a^7/cos(t)^6*sin(t)+5/42*a^7/cos(t)^4*sin(t)+...
%             5/28*a^7/cos(t)^2*sin(t)+5/28*a^7*log(sec(t)+...
%             tan(t))-1/2*a^6/cos(t)^5*sin(t)-2/3*a^6/cos(t)^3*sin(t)-...
%             4/3*a^6/cos(t)*sin(t)+a^5/cos(t)^4*sin(t)+...
%             3/2*a^5/cos(t)^2*sin(t)+3/2*a^5*log(sec(t)+tan(t))-...
%             5/6*a^4/cos(t)^3*sin(t)-5/3*a^4/cos(t)*sin(t)+...
%             1/2*a^2/cos(t)*sin(t);
%         I=W2_int(r0,theta)-W2_int(r0,0);
%
% Note that "W2_int(r0,0)=0".

a=r0; t=theta;

I=2/21*a^7/cos(t)^6*sin(t)+5/42*a^7/cos(t)^4*sin(t)+...
    5/28*a^7/cos(t)^2*sin(t)+5/28*a^7*log(sec(t)+...
    tan(t))-1/2*a^6/cos(t)^5*sin(t)-2/3*a^6/cos(t)^3*sin(t)-...
    4/3*a^6/cos(t)*sin(t)+a^5/cos(t)^4*sin(t)+...
    3/2*a^5/cos(t)^2*sin(t)+3/2*a^5*log(sec(t)+tan(t))-...
    5/6*a^4/cos(t)^3*sin(t)-5/3*a^4/cos(t)*sin(t)+...
    1/2*a^2/cos(t)*sin(t);

I=RBF_scale^2*I;











function I=TPS_integral(misOA,misHA,misHO,RBF_scale)

%--------------------------------------------------------------------------
% Object.
%--------------------------------------------------------------------------
% It computes the integral of the TPS RBF with center the origin,
% RBF scale equal to "RBF_scale", on the triangle OHA, right in H.
% The TPS function is phi(r)=r^2*log(r).
%--------------------------------------------------------------------------

method=2;

switch method
    
    
    
    case 1 % method from cartesian coordinates (working)
        
        % ....................... some assignments ........................
        
        a=misOA; b=misHA; c=misHO; sig=RBF_scale;
        bs=b/sig; cs=c/sig; gam=c/b;
        
        % ...............................................................
        % Note:
        % The primitive is defined as the sum of two functions, i.e. "fI"
        % and "fII", that must be evaluated in "bs" and "0".
        % ...............................................................
        
        % .................... first function analysis ....................
        
        % The code below actually does:
        %
        % fI=@(s,t) -1/18*(s^2+t^2)*(s*t*(5-3*log(s^2+t^2+
        %             eps*(1-abs(sign(t)))))+...
        %     3*(s^2-t^2)*atan(s/(t+eps*(1-sign(s^2+t^2)))));
        % term(1)=(fI(cs,bs)-fI(cs,0));
        
        s=cs; t=bs;
        t2=-1/18*(s^2+t^2)*(s*t*(5-3*log(s^2+t^2+eps*(1-abs(sign(t)))))+...
            3*(s^2-t^2)*atan(s/(t+eps*(1-sign(s^2+t^2)))));
        
        s=cs; t=0;
        t1=-1/18*(s^2+t^2)*(s*t*(5-3*log(s^2+t^2+eps*(1-abs(sign(t)))))+...
            3*(s^2-t^2)*atan(s/(t+eps*(1-sign(s^2+t^2)))));
        
        term(1)=t2-t1;
        
        % .................... second function analysis ....................
        
        % The code below actually does:
        %
        % fII=@(s,t) (1/144)*t^4*(s*(6*(s^2+3)*...
        %     log((s^2+1)*t^2+eps*(1-abs(sign(t))))-7*s^2-33)+24*atan(s));
        % term(2)=(fII(gam,bs)-fII(gam,0));
        
        s=gam; t=bs;
        t2=(1/144)*t^4*(s*(6*(s^2+3)*...
            log((s^2+1)*t^2+eps*(1-abs(sign(t))))-7*s^2-33)+24*atan(s));
        
        s=gam; t=0;
        t1=(1/144)*t^4*(s*(6*(s^2+3)*...
            log((s^2+1)*t^2+eps*(1-abs(sign(t))))-7*s^2-33)+24*atan(s));
        
        term(2)=t2-t1;
        
        % .................... integral on triangle OHA ...................
        
        I=sig^2*(term(1)-term(2));
        
        
        
    case 2 % method from polar coordinates (working)
        
        r0=misHO/RBF_scale; r1=misOA/RBF_scale; theta=acos(r0/r1);
        
        % the code below performs (more quickly!):
        %
        % f=@(t,c) (1/144)*c^4*(12*(cos(2*t) + 2)*...
        %           tan(t)*(sec(t))^2*log(c*sec(t)) +24*t - 26*tan(t)+ ...
        %             (-7)*tan(t)*(sec(t))^2);
        % I=f(r0,theta)-f(r0,0);
        %
        % Note that "f(r0,0)=0" theoretically, but in the above evaluation
        % gives "NaN". By a test in a previous routine, "theta > 0" and
        % "r0>0".
        
        c=r0; t=theta;
        
        I=(1/144)*c^4*(12*(cos(2*t) + 2)*...
            tan(t)*(sec(t))^2*log(c*sec(t)) +24*t - 26*tan(t)+ ...
            (-7)*tan(t)*(sec(t))^2);
        
        I=RBF_scale^2*I;
        
        
        
    case 3 % method from polar coordinates (2006)
        
        r0=misHO;
        r1=misOA;
        
        theta=acos(r0./r1);
        t=tan(theta/2);
        delta=RBF_scale;
        a1=r0./delta;
        t2=t.^2;
        t3=t.^3;
        t4=t.^4;
        t5=t.^5;
        t6=t.^6;
        term1=(a1.^4)./ ((1-t2).^3);
        term2=t.*(3-2*t2+3*t4);
        term3=log(a1.*(1+t2)./(1-t2));
        term4=33*t-38*t3+33*t5;
        term5=-2*(1-3*t2+3*t4-t6).*atan(1./t);
        y2=(1/6)*term1.*(term2.*term3-(1/12)*term4+term5);
        
        term6=(a1./cos(theta)).^3;
        term7=a1./((1+t2).^3);
        Ftheta=(1/6)*term6.*term7.*(term2.*term3-(1/12)*term4+term5);
        F0=-(1/6).*(a1.^4)*pi;
        y=Ftheta-F0;
        
        I=RBF_scale^2*y;
        
end











function I=R3_integral(misOA,misHA,misHO,RBF_scale)

%--------------------------------------------------------------------------
% Object.
%--------------------------------------------------------------------------
% It computes the integral of the "r^3" RBF with center the origin "O",
% RBF scale equal to "RBF_scale", on the triangle OHA, right in H.
%--------------------------------------------------------------------------

method=2;

switch method
    
    case 1 % method from cartesian coordinates (slower)
        
        % ....................... some assignments ........................
        
        a=misOA; b=misHA; c=misHO; sig=RBF_scale;
        bs=b/sig; cs=c/sig; gam=c/b;
        
        % ................................................................
        % Note:
        % The primitive is defined as the sum of two functions, i.e. "fI"
        % and "fII", that must be evaluated in "bs" and "0".
        % ................................................................
        
        % .................... first function analysis ....................
        
        fI=@(s,t) (1/200)*(...
            15*t^5*log(sqrt(s^2+t^2)+s)+...
            35*s*t^3*sqrt(s^2+t^2)+...
            15*s^5*log(sqrt(s^2+t^2)+t)+...
            35*s^3*t*sqrt(s^2+t^2)+...
            (-3)*t^5 );
        
        term(1)=(fI(cs,bs)-fI(cs,0));
        
        % .................... second function analysis ...................
        
        fII=@(s,t) (1/200)*t^4*(...
            5*s*(2*s^2+5)*sqrt((s^2+1)*t^2)+...
            15*t*log(sqrt((s^2+1)*t^2)+s*t+eps*(1-abs(sign(t))))+...
            (-3)*t );
        
        term(2)=(fII(gam,bs)-fII(gam,0));
        
        % .................... integral on triangle OHA ...................
        
        I=sig^2*(term(1)-term(2));
        
    case 2 % method from polar coordinates (faster)
        
        r0=misHO/RBF_scale; r1=misOA/RBF_scale; theta=acos(r0/r1);
        
        % the code below performs (more quickly!):
        %
        % f=@(t,c) (1/80)*c^5*((1/2)*(11*sin(t) + ...
        %      3*sin(3*t))* (sec(t))^4 - 6*log(cos(t/2) - sin(t/2)) + ...
        %      6*log(sin(t/2) + cos(t/2)));
        % I=f(r0,theta)-f(r0,0);
        %
        % Note that "f(r0,0)=0".
        
        c=r0; t=theta;
        
        I=(1/80)*c^5*((1/2)*(11*sin(t) + ...
            3*sin(3*t))* (sec(t))^4 - 6*log(cos(t/2) - sin(t/2)) + ...
            6*log(sin(t/2) + cos(t/2)));
        
        I=RBF_scale^2*I;
        
end











function I=R5_integral(misOA,misHA,misHO,RBF_scale)

%--------------------------------------------------------------------------
% Object.
%--------------------------------------------------------------------------
% It computes the integral of the "phi(r)=r^5" RBF with center in O,
% RBF scale equal to "RBF_scale", on the triangle OHA, right in H.
%--------------------------------------------------------------------------

method=2;

switch method
    case 1 % method from cartesian coordinates (slower, but working)
        
        % ....................... some assignments ........................
        
        a=misOA; b=misHA; c=misHO; sig=RBF_scale;
        bs=b/sig; cs=c/sig; gam=c/b;
        
        % .................................................................
        % Note:
        % The primitive is defined as the sum of two functions, i.e. "fI"
        % and "fII", that must be evaluated in "bs" and "0".
        % .................................................................
        
        % .................... first function analysis ....................
        
        fI=@(s,t) (105*t^7*log(sqrt(s^2 + t^2) + s) + ...
            287*s*t^5*sqrt(s^2 + t^2)+105*s^7*log(sqrt(s^2 + t^2) + t) +...
            287*s^5*t*sqrt(s^2 + t^2)+364*s^3*t^3*sqrt(s^2 + t^2) -...
            15*t^7)/2352;
        
        term(1)=(fI(cs,bs)-fI(cs,0));
        
        % .................... second function analysis ...................
        
        fII=@(s,t) (t^6*(105*t*log(sqrt((s^2 + 1)*t^2) + ...
            +s*t+eps*(1-abs(sign(t)))) + ...
            7*s*(8*s^4 + 26*s^2 + 33) * sqrt((s^2 + 1)* t^2) - 15*t))/2352;
        
        term(2)=(fII(gam,bs)-fII(gam,0));
        
        % .................... integral on triangle OHA ...................
        
        I=sig^2*(term(1)-term(2));
        
        
    case 2 % method from polar coordinates (faster and working)
        
        r0=misHO/RBF_scale; r1=misOA/RBF_scale; theta=acos(r0/r1);
        
        % the code below performs (more quickly!):
        %
        % f=@(t,c) (1/672)*c^7*((1/8)*(198*sin(t) + 85*sin(3*t) + ...
        % 15*sin(5*t))*(sec(t))^6 - 30*log(cos(t/2) - sin(t/2)) + ...
        % 30*log(sin(t/2) + cos(t/2)));
        % I=f(theta,r0)-f(0,r0);
        %
        % Note that "f(0,r0)=0".
        
        c=r0; t=theta;
        
        I=(1/672)*c^7*((1/8)*(198*sin(t) + 85*sin(3*t) + ...
            15*sin(5*t))*(sec(t))^6 - 30*log(cos(t/2) - sin(t/2)) + ...
            30*log(sin(t/2) + cos(t/2)));
        
        I=RBF_scale^2*I;
        
end











function I=R7_integral(misOA,misHA,misHO,RBF_scale)

%--------------------------------------------------------------------------
% Object.
%--------------------------------------------------------------------------
% It computes the integral of the "phi(r)=r^5" RBF with center in O,
% RBF scale equal to "RBF_scale", on the triangle OHA, right in H.
%--------------------------------------------------------------------------

method=2;

switch method
    case 1 % method from cartesian coordinates (slower, but working)
        error('\n \t method not implemented');
        
        
    case 2 % method from polar coordinates (faster and working)
        
        r0=misHO/RBF_scale; r1=misOA/RBF_scale; theta=acos(r0/r1);
        
        % the code below performs (more quickly!):
        %
        % f=@(t,c) (c^9* (10106*tan(t)*(sec(t))^7 - 13440*...
        %   log(cos(t/2) - sin(t/2)) + ...
        %   7*(sec(t))^8* (766* sin(3* t) + 230* sin(5* t) + ...
        %   30* sin(7* t) +840*cos(2* t)* log(sin(t/2) + cos(t/2)) +  ...
        %   420* cos(4* t)* log(sin(t/2) + cos(t/2)) + 120* cos(6* t)* ...
        %   log(sin(t/2) + cos(t/2)) + 15* cos(8* t)* log(sin(t/2) + ...
        %   cos(t/2)) +525* log(sin(t/2) + cos(t/2)))))/442368;
        % I=f(theta,r0)-f(0,r0);
        %
        % Note that "f(0,r0)=0".
        
        c=r0; t=theta;
        
        I=(c^9* (10106*tan(t)*(sec(t))^7 - 13440*log(cos(t/2) - sin(t/2)) + ...
            7*(sec(t))^8* (766* sin(3* t) + 230* sin(5* t) + 30* sin(7* t) + ...
            840*cos(2* t)* log(sin(t/2) + cos(t/2)) + 420* cos(4* t)* ...
            log(sin(t/2) + cos(t/2)) + 120* cos(6* t)* log(sin(t/2) + cos(t/2)) + ...
            15* cos(8* t)* log(sin(t/2) + cos(t/2)) + ...
            525* log(sin(t/2) + cos(t/2)))))/442368;
        
        
        I=RBF_scale^2*I;
        
end













function I=W_integral_19(misOA,misHA,misHO,RBF_scale,k)

%--------------------------------------------------------------------------
% Object.
%--------------------------------------------------------------------------
% It computes the integral of the Wendland RBF (W8)
%         phi=@(r) ((1-r))^10.*(429*r^4 + 450*r^3 + 210*r^2 + 5)
% with center in O, RBF scale equal to "RBF_scale", on the
% triangle OHA, right in H.
%--------------------------------------------------------------------------

% ....................... some assignments ................................

a=misOA; b=misHA; c=misHO; sig=RBF_scale;
bs=b/sig; cs=c/sig; gam=c/b;

% measure in R^2 of the Wendland RBF with scale equal to 1.

switch k
    case 0
        mu0=pi/6;
    case 4
        mu0=pi/3;
    case 6
        mu0=14*pi/156;
    case 8
        mu0=(3*pi)/8;
end
% ..........  the region is union of a sector and a right triangle ........

if RBF_scale < a & RBF_scale > c
    
    % ... computing measure of sector ...
    
    thetaM=asin(c/sig);
    thetam=asin(c/a);
    mu_sector=sig^2*mu0*(thetaM-thetam)/(2*pi);
    
    % ... computing RBF integral on a triangle, in the support ...
    
    % we reduce the computation to a right triangle, OHJ, right in H,
    % with "mis0J=sig".
    misOJ=RBF_scale;
    misHJ=sqrt(RBF_scale^2-c^2);
    int_insupp =W_integral_19_insupp(misOJ,misHJ,misHO,RBF_scale,k);
    
    I=mu_sector+int_insupp;
    return;
    
end

% ..........  the region is a sector ..........

% ... in this instance the radius of the sector is RBF_scale ...
if c >= RBF_scale
    
    thetaM=pi/2;
    thetam=asin(c/a);
    I=sig^2*mu0*(thetaM-thetam)/(2*pi);
    
    
else
    % ..........  the region is a right triangle ..........
    I=W_integral_19_insupp(misOA,misHA,misHO,RBF_scale,k);
end










function I=MW0_integral_19(misOA,misHA,misHO,RBF_scale)

%--------------------------------------------------------------------------
% Object.
%--------------------------------------------------------------------------
% It computes the integral of the missing Wendland RBF (W2) with center the
% origin, RBF scale equal to "RBF_scale", on the triangle OHA, right in H.
%--------------------------------------------------------------------------

% ....................... some assignments ................................
a=misOA; b=misHA; c=misHO; sig=RBF_scale;
bs=b/sig; cs=c/sig; gam=c/b;

mu0=(sqrt(2*pi)/15); % measure in R^2 of the Wendland RBF with scale equal to 1.
% mu0=1.014851768682233e-01;

% ..........  the region is union of a sector and a right triangle ........

if RBF_scale < a & RBF_scale > c
    
    % ... computing measure of sector ...
    
    thetaM=asin(c/sig);
    thetam=asin(c/a);
    mu_sector=sig^2*mu0*(thetaM-thetam)/(2*pi);
    
    % ... computing RBF integral on a triangle, in the support ...
    
    % we reduce the computation to a right triangle, OHJ, right in H,
    % with "mis0J=sig".
    misOJ=RBF_scale;
    misHJ=sqrt(RBF_scale^2-c^2);
    int_insupp=MW0_integral_19_insupp(misOJ,misHJ,misHO,RBF_scale);
    I=mu_sector+int_insupp;
    return;
    
end

% ..........  the region is a sector ..........

% ... in this instance the radius of the sector is RBF_scale ...
if c >= RBF_scale
    
    thetaM=pi/2;
    thetam=asin(c/a);
    I=sig^2*mu0*(thetaM-thetam)/(2*pi);
    
    
else
    % ..........  the region is a right triangle ..........
    I=MW0_integral_19_insupp(misOA,misHA,misHO,RBF_scale);
end





function I=MW0_integral_19_insupp(misOA,misHA,misHO,RBF_scale,method)

% Note:
% The three codes below are more or less equivalent, they provide same
% results, but something is wrong in the procedure ...

r0=misHO/RBF_scale; r1=misOA/RBF_scale; theta=acos(r0/r1);

if nargin<5, method=1; end

switch method
    case 1 % polar, with first primitive
        
        pert=@(r) eps*(1-abs(sign(r)));
        H_primitive=@(r) 3*sqrt(1-r.^2).*r.^2/(10*sqrt(2*pi))-...
            sqrt(1-r.^2)/(15*sqrt(2*pi))+...
            2/15*sqrt(2/pi)*sqrt(1-r.^2).*r.^4+...
            r.^4.*log((r+pert(r))./(sqrt(1-r.^2)+1))/...
            (2*sqrt(2*pi));
        
        integrand=@(t) H_primitive(r0./cos(t));
        
        I1=integral(integrand,0,theta,'AbsTol',10^(-14),'RelTol',10^(-14));
        I0=theta*(-2.659615202676218e-02);
        I=RBF_scale^2*(I1-I0);
        
    case 2 % polar, full numerics
        
        pert=@(r) eps*(1-abs(sign(r)));
        phi=@(r) (sqrt(2)/(3*sqrt(pi)))*...
            (...
            3*r.^2.*log(r./(1+sqrt(1-r.^2))+pert(r))+...
            (2*r.^2+1).*sqrt(1-r.^2)...
            );
        
        phir=@(t,r) r.*phi(r)+0*t;
        espo=@(s) r0./cos(s);
        Inum=integral2(phir,0,theta,0,espo,'AbsTol',10^(-15));
        I=RBF_scale^2*Inum;
        
    case 3 % cartesian, numerics
        
        sqrtxy=@(x,y) sqrt(x.^2+y.^2);
        pert=@(r) eps*(1-abs(sign(r)));
        phi=@(r) sqrt(2)/(3*sqrt(pi)) *...
            (...
            3*r.^2.*log(r./(1+sqrt(1-r.^2))+pert(r))+...
            (2*r.^2+1).*sqrt(1-r.^2) ...
            );
        phi_cart=@(x,y) phi(sqrtxy(x,y));
        
        b=sqrt(r1^2-r0^2);
        ymin=@(x) (r0/b)*x;
        
        Inum=integral2(phi_cart,0,b/RBF_scale,ymin,r0/RBF_scale,...
            'AbsTol',10^(-15));
        I=RBF_scale^2*Inum;
        
end












function I_tri=W_integral_19_insupp(misOA,misHA,misHO,RBF_scale,k)

r0=misHO/RBF_scale; r1=misOA/RBF_scale; theta=acos(r0/r1);
c=r0; t=theta;

% These coefficients are obtained by symbolic calculus, integrating
% "phi(r)*r" in Maple.
% Example (Wendland W0):
% syms r
% int(r*(1-r)^2)
% ans =
% (r^2*(3*r^2 - 8*r + 6))/12
% expand((r^2*(3*r^2 - 8*r + 6))/12)
% ans =
% r^4/4 - (2*r^3)/3 + r^2/2

switch k
    case 0
        coeffs=[0 0 1/2 -2/3 1/4];
    case 2
        coeffs=[0 0 1/2 0 -5/2 4 -5/2 4/7];
    case 4
        coeffs=[0 0 3/2 0 -7 0 35 -64 105/2 -64/3 7/2];
    case 6
        coeffs=[0 0 1/2 0 -11/4 0 11 0 -231/4 352/3 -231/2 64 -77/4 32/13];
    case 8
        coeffs=[429/16 -256 2145/2 -2560 15015/4 -3328 3003/2 ...
            0 -2145/8 0 143/2 0 -65/4 0 5/2 0 0];
        coeffs=fliplr(coeffs);
        
    otherwise
        error('case not implemented');
end

L=length(coeffs);
coeffs=coeffs.*(c.^(0:L-1));

I=[t; log(abs(tan((t/2)+(pi/4))))];

for j=3:L
    n=j-1; % degree of (1/cos(t))^n to be computed
    I(j,1)=sin(t)/((n-1)*(cos(t))^(n-1))+((n-2)/(n-1))*I(j-2);
end

I_tri=RBF_scale^2*(coeffs*I);






function I=M1_integral_19(misOA,misHA,misHO,RBF_scale)

% Note:
% The code below resorts to 1D adaptive routine, to compute the required
% bivariate integral.
% The RBF is the Matern "phi=@(r) exp(-r)" obtained for "beta_1=3/2".

r0=misHO/RBF_scale; r1=misOA/RBF_scale; theta=acos(r0/r1);
psi=@(t) (-exp(-t).*(t + 1))+1;
psi_eps=@(t) psi(r0./(RBF_scale*cos(t)));
tol=10^(-15);
Ipsi=integral(psi_eps,0,theta,'AbsTol',tol,'RelTol',tol);
I=RBF_scale^2*Ipsi;






function I=M2_integral_19(misOA,misHA,misHO,RBF_scale)

% Note:
% The code below resorts to 1D adaptive routine, to compute the required
% bivariate integral.
% The RBF is the Matern "phi=@(r) (1+r)*exp(-r)" obtained for "beta_2=5/2".

r0=misHO/RBF_scale; r1=misOA/RBF_scale; theta=acos(r0/r1);
psi=@(r) (-exp(-r).*(r.^2 + 3*r + 3))+3;
psi_eps=@(t) psi(r0./(RBF_scale*cos(t)));
tol=10^(-15);
Ipsi=integral(psi_eps,0,theta,'AbsTol',tol,'RelTol',tol);
I=RBF_scale^2*Ipsi;







function pgon=polygon_to_polyshape(xv,yv,iv)

S=class(xv);

S_double = strcmp(S,'double');

if S_double % distinguishing polyshape class from double.
    if nargin < 2
        pts=xv;
        if size(pts,1) <= size(pts,2) % working with column vectors.
            pts=xv';
        end
        if pts(1,:) == pts(end,:)
            pts=pts(1:end-1,:);
        end
        xv=pts(:,1); yv=pts(:,2); iv=length(xv);
    else
        if size(xv,1) <= size(xv,2) % working with column vectors.
            xv=xv';
        end
        if size(yv,1) <= size(yv,2) % working with column vectors.
            yv=yv';
        end
        %     if [xv(1) yv(1)] == [xv(end) yv(end)]
        %         xv=xv(1:end-1); yv=yv(1:end-1);
        %     end
        if nargin < 3
            iv=length(xv);
        else
            if isempty(iv)
                iv=length(xv);
            end
        end
    end
    
    % Here we convert the given information into cells of coordinates, where
    % each component represent a geometrical component of the polygon, not
    % connected or at most tangent to other geometrical components.
    
    xvc={xv(1:iv(1))}; yvc={yv(1:iv(1))};
    ii1=1; ii2=iv(1);
    for ii=2:length(iv)
        ii1=ii2+1;
        ii2=ii2+iv(ii);
        xvl=xv(ii1:ii2);
        yvl=yv(ii1:ii2);
        
        xvc{end+1}=xvl; yvc{end+1}=yvl;
    end
    
    pgon = polyshape(xvc,yvc);
else
    pgon=xv;
    [xvc,yvc]=boundary(pgon);
end





