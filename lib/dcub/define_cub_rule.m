
function [XW,dbox,domain_structure]=define_cub_rule(domain_structure,...
    deg,flag_compression)

%--------------------------------------------------------------------------
% Object:
% Define algebraic cubature rule of algebraic degree of precision "deg" on
% domain defined by "domain_structure".
%--------------------------------------------------------------------------
% Input:
% domain_structure: structure defining the domain.
% deg: algebraic degree of precision of the rule.
% flag_compression: if "1" compression will be attempted, otherwise the
%      code will propose an alternative rule.
%--------------------------------------------------------------------------
% XW: cubature rule of degree "deg".
% dbox: bounding box as [xLimit; yLimit] where
%     xLimit=[xmin xmax], yLimit=[ymin ymax].
% domain_structure: domain structure including "domain_structure.dbox" and
%      "domain_structure.area".
%--------------------------------------------------------------------------
% Example of "domain_structure":
%
% Polygon:
%      domain_struct.domain='polygon';
%      domain_struct.vertices=[0.1 0; 0.7 0.2; 1 0.5; 0.75 0.85];
% It is a polygon with vertices described by rows of the value of
% "domain_struct.vertices". The polygon is described counterclockwise by
% vertices.
%
% Lune:
%      domain_struct.domain='lune';
%      x1=0; y1=0; r1=1.5; x2=1; y2=0; r1=1.5;
%      domain_struct.parms=[x1 y1 r1 x2 y2 r2];
%
% Circular Annular Sector:
%     domain_struct.domain='circular-annular-sector'
%     x1=1; x2=1.5; r1=0.5; r2=2;
%     domain_struct.parms=[x1 x2 r1 r2];
%
% Disk:
%     domain_struct.domain='disk'
%         x1=1; x2=1.5; r1=0.5;
%         domain_struct.parms=[x1 x2 r1];
%
% Asymmetric circular sector:
%     domain_struct.domain='asymmetric-circular-sector'
%         a=1; b=1; c=2; d=2; r=3;
%         domain_struct.parms=[a b c d r];
%
% Asymmetric annulus:
%     domain_struct.domain='asymmetric-annulus'
%         a=1; b=1; c=2; d=2; r1=3; r2=4;
%         domain_struct.parms=[a b c d r1 r2];
%
% Vertical Circular zone:
%     domain_struct.domain='vertical-circular-zone'
%         a=1; b=1; r=4; alpha=-pi/4;
%         domain_struct.parms=[a b r alpha];
%
% Horizontal Circular zone:
%     domain_struct.domain='horizontal-circular-zone'
%         a=1; b=1; r=4; alpha=-pi/4;
%         domain_struct.parms=[a b r alpha];
%
% Circular segment:
%     domain_struct.domain='circular-segment'
%         a=1; b=1; r=4; alpha=-pi/4; beta=pi/6;
%         domain_struct.parms=[a b r alpha beta];
%
% Symmetric lens:
%     domain_struct.domain='symmetric-lens'
%         % in this case "domain_parms" is a vector [a r] (with a < r)
%         a=1; r=3;
%         domain_struct.parms=[a r];
%
% Butterfly:
%     domain_struct.domain='butterfly'
%         a=1; b=2; r=3; alpha=-pi/4; beta=pi/6;
%         domain_struct.parms=[a b r alpha beta];
%
% Candy:
%     domain_struct.domain='candy'
%         % in this case "domain_parms" is a vector [a r alpha] with
%         % "-alpha > acos(a/r)"
%         a=0.5; r=1; alpha=-1.5;
%         domain_struct.parms=[a r alpha];
%
%--------------------------------------------------------------------------

% ......................... troubleshooting ...............................

if ~isempty(domain_structure)

    % ................... check input ...................

    if ~isa(domain_structure,'struct')
        error('MATLAB: "define_cub_rule"',...
            getString(message(...
            'Variable "domain_structure" not defined')));
    end

    % ................... define domain ...................

    if isfield(domain_structure,'domain')
        domain = domain_structure.domain;
    else
        error('MATLAB: "define_cub_rule"',...
            getString(message(...
            'Variable "domain" not defined')));
    end

    % ..... domain parms check ......

    S(1)=strcmp(domain,'polygon');
    S(2)=strcmp(domain,'lune');
    S(3)=strcmp(domain,'circular-annular-sector');
    S(4)=strcmp(domain,'disk');
    S(5)=strcmp(domain,'asymmetric-circular-sector');
    S(6)=strcmp(domain,'asymmetric-annulus');
    S(7)=strcmp(domain,'vertical-circular-zone');
    S(8)=strcmp(domain,'horizontal-circular-zone');
    S(9)=strcmp(domain,'circular-segment');
    S(10)=strcmp(domain,'symmetric-lens');
    S(11)=strcmp(domain,'butterfly');
    S(12)=strcmp(domain,'candy');
    S(13)=strcmp(domain,'convex-level-curve');
    S(14)=strcmp(domain,'circle-arc');
    S(15)=strcmp(domain,'sphere');
    S(16)=strcmp(domain,'spherical-rectangle');
    S(17)=strcmp(domain,'pyramid');
    S(18)=strcmp(domain,'3D_rect');
    S(19)=strcmp(domain,'rotation-domain');
    S(20)=strcmp(domain,'polygonal-boundary');
    S(21)=strcmp(domain,'pluriinterval');
    S(22)=strcmp(domain,'sector');
    S(23)=strcmp(domain,'pyramid-boundary');
    S(24)=strcmp(domain,'asymmetric-circular-sector-boundary');
    S(25)=strcmp(domain,'asymmetric-annulus-boundary');
    S(26)=strcmp(domain,'vertical-circular-zone-boundary');
    S(27)=strcmp(domain,'horizontal-circular-zone-boundary');
    S(28)=strcmp(domain,'circular-segment-boundary');
    S(29)=strcmp(domain,'symmetric-lens-boundary');
    S(30)=strcmp(domain,'butterfly-boundary');
    S(31)=strcmp(domain,'candy-boundary');
    S(32)=strcmp(domain,'sector-boundary');
    S(33)=strcmp(domain,'rotation-surface');
    S(34)=strcmp(domain,'spherical-triangle');
    S(35)=strcmp(domain,'NURBS');
    S(36)=strcmp(domain,'union-disks');
    S(37)=strcmp(domain,'square');
    S(38)=strcmp(domain,'unit-square[0,1]x[0,1]');
    S(38)=strcmp(domain,'rectangle');
    S(39)=strcmp(domain,'triangle');
    S(40)=strcmp(domain,'polygcirc');
    S(41)=strcmp(domain,'unit-simplex');
    S(42)=strcmp(domain,'spherical-polygon');

end


if nargin < 2, deg=10; end
if nargin < 3, flag_compression=0; end


% .......................... WAM definition ...............................
dbox=[];

switch domain

    case 'polygon'

        PSHtest=isfield(domain_structure,'polyshape');

        if PSHtest
            do_polyshape=1;
            XY=domain_structure.polyshape;
        else
            do_polyshape=0;
            XY=domain_structure.vertices;
        end

        deg_max_comp=15;

        if do_polyshape
            % ... XY is a polyshape object ...
            if (deg <= deg_max_comp) | (flag_compression == 1)
                [XWF,~,~,~,~,XW]=cub_polygon(deg,XY);
            else
                XW=cub_polygon(deg,XY);
            end
        else
            % ... XY is a matrix ...
            if (deg <= deg_max_comp) | (flag_compression == 1)
                [XWF,~,~,~,~,XW]=cub_polygon(deg,XY(:,1),XY(:,2));
            else
                XW=cub_polygon(deg,XY(:,1),XY(:,2));
            end

        end


    case 'polygonal-boundary'
        error('Not implemented yet.');

        XW=[];


    case 'lune'

        center1=domain_structure.centers(1,:);
        center2=domain_structure.centers(2,:);

        r1=domain_structure.radii(1);
        r2=domain_structure.radii(2);

        x1=center1(1); y1=center1(2);
        x2=center2(1); y2=center2(2);

        XW= cub_lune(deg,x1,y1,r1,x2,y2,r2);

    case 'circular-annular-sector'

        % ...................... decode variables .........................

        a=domain_structure.center(1);
        b=domain_structure.center(2);
        r1=domain_structure.radii(1);
        r2=domain_structure.radii(2);
        alpha=domain_structure.angles(1);
        beta=domain_structure.angles(2);

        % ........................ compute rule ...........................

        if r1 >= r2
            error('MATLAB: "define_cub_rule" on', ...
                'circular-annular-sector',getString(message(...
                'r1 must be smaller than r2')));
        end

        A=[r2 0;r1 0]; B=[0 r2;0 r1]; C=[a b;a b];

        XW= cub_ellblend(deg,A,B,C,alpha,beta);


    case 'disk'

        % ...................... decode variables .........................

        a=domain_structure.center(1);
        b=domain_structure.center(2);
        r=domain_structure.radius;

        % ........................ compute rule ...........................

        XW=cub_disk(deg,[a b],r);

    case 'sector'

        % ...................... decode variables .........................

        a=domain_structure.center(1);
        b=domain_structure.center(2);
        r2=domain_structure.radius; r1=0;
        alpha=domain_structure.angles(1);
        beta=domain_structure.angles(2);

        % ........................ compute rule ...........................

        method=1;

        switch method
            case 1 % cub_ellblend

                A=[r2 0;r1 0]; B=[0 r2;0 r1]; C=[a b;a b];
                XW = cub_ellblend(deg,A,B,C,alpha,beta);

            case 2 % cub_circsect

                if beta < alpha, beta=beta+2*pi; end

                omega=(beta-alpha)/2;
                XW0=cub_circsect(deg,omega,r1,r2);
                th=(beta+alpha)/2;
                rotmat=[cos(th) -sin(th); sin(th) cos(th)];
                X1=(rotmat*(XW0(:,1:2))')';
                XW=[X1(:,1)+a X1(:,2)+b XW0(:,3)];

        end

    case 'asymmetric-circular-sector'

        % ...................... decode variables .........................

        centers=domain_structure.centers;
        r=domain_structure.radius;
        anglesV=domain_structure.angles;

        a=centers(1,1); b=centers(1,2);
        c=centers(2,1); d=centers(2,2);

        alpha=anglesV(1); beta=anglesV(2);

        % ........................ compute rule ...........................

        A=[r 0;0 0]; B=[0 r;0 0]; C=[a b;c d];
        XW= cub_ellblend(deg,A,B,C,alpha,beta);


    case 'asymmetric-annulus'

        % ...................... decode variables .........................

        centers=domain_structure.centers;
        rV=domain_structure.radii;

        a=centers(1,1); b=centers(1,2);
        c=centers(2,1); d=centers(2,2);

        r1=rV(1); r2=rV(2);

        % ........................ compute rule ...........................

        A=[r1 0; r2 0]; B=[0 r1;0 r2]; C=[a b;c d];
        alpha=-pi; beta=pi;
        XW= cub_ellblend(deg,A,B,C,alpha,beta);


    case 'vertical-circular-zone'

        % ...................... decode variables .........................

        a=domain_structure.center(1);
        b=domain_structure.center(2);
        r=domain_structure.radius;
        alpha=domain_structure.angles(1);
        beta=domain_structure.angles(2);

        % ........................ compute rule ...........................

        A=[r 0;r 0]; B=[0 r;0 -r]; C=[a b;a b];
        XW= cub_ellblend(deg,A,B,C,alpha,beta);


    case 'horizontal-circular-zone'

        % ...................... decode variables .........................

        a=domain_structure.center(1);
        b=domain_structure.center(2);
        r=domain_structure.radius;
        alpha=domain_structure.angles(1);
        beta=domain_structure.angles(2);

        % ........................ compute rule ...........................

        A=[r 0;-r 0]; B=[0 r;0 r]; C=[a b;a b];

        XW= cub_ellblend(deg,A,B,C,alpha,beta);


    case 'circular-segment'

        % ...................... decode variables .........................

        a=domain_structure.center(1);
        b=domain_structure.center(2);
        r=domain_structure.radius;
        alpha=domain_structure.angles(1);
        beta=domain_structure.angles(2);

        % ........................ compute rule ...........................

        if beta < alpha, beta=beta+2*pi; end

        omega=(beta-alpha)/2;

        XW0 = cub_circsegm(deg,omega,r);
        [THref,R]=cart2pol(XW0(:,1),XW0(:,2));

        TH=(beta+alpha)/2+THref;

        XW=[a+R.*cos(TH) b+R.*sin(TH) XW0(:,3)];


    case 'symmetric-lens'

        % ...................... decode variables .........................

        a=domain_structure.center(1);
        r=domain_structure.radius;

        % ........................ compute rule ...........................

        if a > r
            error('MATLAB: "define_cub_rule"',getString(message(...
                'a must be strictly smaller than b')));
        end

        beta=acos(a/r); alpha=-beta;
        A=[r 0;-r 0]; B=[0 r;0 r]; C=[-a 0;a 0];
        XW= cub_ellblend(deg,A,B,C,alpha,beta);


    case 'butterfly'

        % ...................... decode variables .........................

        a=domain_structure.center(1);
        b=domain_structure.center(2);
        r=domain_structure.radius;
        alpha=domain_structure.angles(1);
        beta=domain_structure.angles(2);

        % ........................ compute rule ...........................

        A=[r 0;-r 0]; B=[0 r;0 -r]; C=[a b;a b];
        XW= cub_ellblend(deg,A,B,C,alpha,beta);



    case 'candy'

        % ...................... decode variables .........................

        a=domain_structure.center(1);
        r=domain_structure.radius;
        alpha=domain_structure.angle(1);
        beta=-alpha;

        % ........................ compute rule ...........................

        acos_a_r=acos(a/r);
        if beta <= acos_a_r
            error('MATLAB: "define_cub_rule"',getString(message(...
                'it must be "-alpha > acos(a/r)"')));
        end
        
        A=[r 0;-r 0]; B=[0 r;0 r]; C=[-a 0;a 0];
        XW=cub_ellblend(deg,A,B,C,alpha,beta);

    case 'sphere'

        XW=cub_sphere(deg);

    case 'spherical-rectangle'

        sph_rect_parms=domain_structure.angles; % [0 pi; 0 2*pi];
        [~,XW]=cub_sphrect(deg,sph_rect_parms);

    case 'spherical-triangle'

        sphtri_vertices=domain_structure.vertices;
        XW = cub_sphtri(deg,sphtri_vertices(1,:),sphtri_vertices(2,:),...
            sphtri_vertices(3,:));

    case 'spherical-polygon'

        sphpgon_vertices=domain_structure.vertices;
        XW=cub_sphpgon(deg,sphpgon_vertices);
        xLimit=[min(XW(:,1)) max(XW(:,1))];
        yLimit=[min(XW(:,2)) max(XW(:,2))];
        zLimit=[min(XW(:,3)) max(XW(:,3))];

        dbox=[xLimit; yLimit; zLimit];

    case 'NURBS'

        [XW,res] = cub_NURBS(deg,domain_structure.NURBS);
        
        W=XW(:,3);
        ineg=find(W < 0);

        if isempty(ineg) == 0
            fprintf('\n \t * Negative weights: %5.0f over %5.0f',...
                length(ineg),length(W));
        end

        if res > 10^(-12)
            fprintf('\n \t * Moment Residual: 1.3e',res);
        end


    case 'union-disks'

        % ...................... decode variables .........................

        centerV=domain_structure.centers;
        rV=domain_structure.radii;

        % ........................ compute rule ...........................

        [~,XW]=cub_uniondisks(centerV,rV,deg,1);

    case 'square'

        % ........................ compute rule ...........................

        XW=cub_square(deg);

    case {'unit-square[0,1]x[0,1]','rectangle'}

        % ........................ compute rule ...........................

        % rule in [-1,1]^2

        XWr=cub_square(deg);
        bbox=domain_structure.dbox;
        xLimit=bbox(1,:);
        yLimit=bbox(2,:);

        % scaling rule
        Xr=XWr(:,1); a=xLimit(1); b=xLimit(2);
        X=(a+b)/2+(b-a)*Xr/2;

        Yr=XWr(:,2); c=yLimit(1); d=yLimit(2);
        Y=(d+c)/2+(d-c)*Yr/2;

        Wr=XWr(:,3);
        area_rect=diff(xLimit)*diff(yLimit);

        W=(area_rect/4)*Wr;

        XW=[X Y W];

    case 'triangle'

        % ........................ compute rule ...........................

        Pgon=domain_structure.polyshape;

        vertices=Pgon.Vertices;

        xw=cub_triangle(deg);
        x=xw(:,1); y=xw(:,2); ww=xw(:,3);
        pts_bar=[x y 1-x-y];

        xx=pts_bar*vertices(1:3,:);
        ww=(area(Pgon)/0.5)*ww;

        XW=[xx ww];

    case 'polygcirc'

        % ...................... decode variables .........................

        ab=domain_structure.extrema;
        a=ab(1,:); b=ab(2,:);
        v=domain_structure.vertices;
        cc=domain_structure.center;
        r=domain_structure.radius;
        conv=domain_structure.convex;

        % ........................ compute rule ...........................

        XW = cub_polygcirc(deg,v,a,b,cc,r,conv);

        L=size(XW,1);
        dim_polyspace=(deg+1)*(deg+2)/2;

        if L > dim_polyspace & (dim_polyspace <= 1000)
            X=XW(:,1:2); u=XW(:,3);
            [XC,uc,momerr] = dCATCH(deg,X,u,[],[],dim_polyspace);
            XW=[XC uc];
        end


    case 'unit-simplex'

        % ........................ compute rule ...........................

        XW=cub_triangle(deg);

end

domain_structure.area=sum(XW(:,3));

if nargout > 1
    if isempty(dbox)
        dbox=domain_boundingbox(domain_structure);
    end
end




