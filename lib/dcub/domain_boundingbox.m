

function domain_structure_dbox=domain_boundingbox(domain_structure)

% domain_structure is a structure containing at least:
% 1. domain_structure.domain: string with domain name
% 2. domain_structure.parms : domain structure (variable from domain).

domain_structure_dbox=[];
domain_string=domain_structure.domain;



switch domain_string

    case {'polygon','square','rectangle','triangle',...
            'unit-square[0,1]x[0,1]'}

        [xlimit,ylimit]=boundingbox(domain_structure.polyshape);
        domain_structure_dbox=[xlimit; ylimit];


    case {'lune','circular-annular-sector','sector',...
            'vertical-circular-zone',...
            'horizontal-circular-zone','circular-segment','butterfly'}

        domain_structure_dbox=compute_bbox_NURBS(domain_structure);

    case 'disk'

        % ...................... decode variables .........................
        center=domain_structure.center;
        r=domain_structure.radius;

        % ................... determine bounding box ......................
        xlimit=[center(1)-r center(1)+ r];
        ylimit=[center(2)-r center(2)+ r];

        domain_structure_dbox=[xlimit; ylimit];

    case 'asymmetric-circular-sector'

        % ...................... decode variables .........................

        structure_RS=domain_structure.NURBS;
        [xyw, ~, ~, ~, ~,bbox] = cub_NURBS(0,structure_RS);

        % ................... determine bounding box ......................

        xLimit=[bbox(1), bbox(2)]; yLimit=[bbox(3), bbox(4)];
        domain_structure_dbox=[xLimit; yLimit];

    case 'asymmetric-annulus'

        % ...................... decode variables .........................

        centerV=domain_structure.centers;
        rV=domain_structure.radii;

        a=centerV(1,1); b=centerV(1,2); c=centerV(2,1); d=centerV(2,2);
        r1=rV(1); r2=rV(2);

        % ................... determine bounding box ......................
        
        xLimit=[min(a-r1,c-r2) max(a+r1,c+r2)];
        yLimit=[min(b-r1,d-r2) max(b+r1,d+r2)];
        domain_structure_dbox=[xLimit; yLimit];


    case 'symmetric-lens'

        % ....................... decode variables .......................
        a=domain_structure.center; a=abs(a);
        r=domain_structure.radius;
        
        center1=[-a 0]; center2=[a 0];

        % ................... determine bounding box ......................
        xlimit=[a-r -a+ r];
        c=sqrt(r^2-a^2); ylimit=[-c c];

        domain_structure_dbox=[xlimit; ylimit];



    case 'candy'

        % ...................... decode variables .........................

        a=domain_structure.center;
        r=domain_structure.radius;
        alpha=domain_structure.angle;

        % ................... determine bounding box ......................

        domain_structure_dbox=[];

    case 'NURBS'

        % ...................... decode variables .........................

        structure_RS=domain_structure.NURBS;
        [xyw, ~, ~, ~, ~,bbox] = cub_NURBS(0,structure_RS);

        % ................... determine bounding box ......................

        xLimit=[bbox(1), bbox(2)]; yLimit=[bbox(3), bbox(4)];
        domain_structure_dbox=[xLimit; yLimit];


    case 'union-disks'

        % ...................... decode variables .........................

        centers=domain_structure.centers;
        rs=domain_structure.radii;

        % ................... determine bounding box ......................
        Xc=centers(:,1); Yc=centers(:,2);

        xLimits=[min(Xc-rs) max(Xc+rs)];
        yLimits=[min(Yc-rs) max(Yc+rs)];
        domain_structure_dbox=[xLimits; yLimits];

    case 'polygcirc'

        % ...................... decode variables .........................

        a=domain_structure.extrema(1,:); 
        b=domain_structure.extrema(2,:);
        v=domain_structure.vertices;
        cc=domain_structure.center;
        r=domain_structure.radius;
        conv=domain_structure.convex;

        % ................... determine bounding box ......................

        vertices=[a; v; b];
        X=[a(1); v(:,1); b(1)]; Y=[a(2); v(:,2); b(2)];

        if conv
            X=[X; cc(1)+r; cc(1)-r]; Y=[Y; cc(2)+r; cc(2)-r];
        end

        xLimits=[min(X) max(X)]; yLimits=[min(Y) max(Y)];
        domain_structure_dbox=[xLimits; yLimits];

    case 'unit-simplex'

        % ................... determine bounding box ......................

        xLimits=[0 1]; yLimits=[0 1];
        domain_structure_dbox=[xLimits; yLimits];

end








function domain_structure_dbox=compute_bbox_NURBS(domain_structure)

[~,~,~,~,~,~,boxVx]=indomain_NURBS([],domain_structure.NURBS);

% rectangle is [bbox(1), bbox(2)] x [bbox(3), bbox(4)];
xLimit=[min(boxVx(:,1)), max(boxVx(:,2))];
yLimit=[min(boxVx(:,3)), max(boxVx(:,4))];

domain_structure_dbox=[xLimit; yLimit];


