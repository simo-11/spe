
function in=indomain_routine(domain_structure,pts)

% 1: 'polygon';
% 2: 'disk';
% 3: 'lune';
% 4: 'circular-annular-sector';
% 5: 'sector';
% 6: 'asymmetric-circular-sector';
% 7: 'asymmetric-annulus';
% 8: 'vertical-circular-zone';
% 9: 'horizontal-circular-zone';
% 10: 'circular-segment';
% 11: 'symmetric-lens';
% 12: 'butterfly';            % NOT PASSED!
% 13: 'candy';                % NOT PASSED!
% 14: 'NURBS';
% 15: 'union-disks';          % 10^(-11)
% 16: 'asymmetric-circular-sector';

domain_str=domain_structure.domain;

switch domain_str

    case {'polygon','unit-simplex','square','rectangle','triangle', ...
            'unit-square[0,1]x[0,1]'}

        % ...................... indomain code ............................

        P=domain_structure.polyshape;
        X=pts(:,1); Y=pts(:,2);
        [TFin,TFon]=isinterior(P,X,Y);
        in=TFin & not(TFon);

    case 'disk'

        % ...................... decode variables .........................

        center=domain_structure.center;
        r=domain_structure.radius;

        % ...................... indomain code ............................

        dist_C_pts=sqrt((pts(:,1)-center(1)).^2+...
            (pts(:,2)-center(2)).^2);
        in=(dist_C_pts <= r);

    case 'lune'

        % ...................... decode variables .........................

        center1=domain_structure.centers(1,:);
        center2=domain_structure.centers(2,:);
        r1=domain_structure.radii(1);
        r2=domain_structure.radii(2);

        dist1=sqrt((pts(:,1)-center1(1)).^2+...
            (pts(:,2)-center1(2)).^2);

        dist2=sqrt((pts(:,1)-center2(1)).^2+...
            (pts(:,2)-center2(2)).^2);

        in=(dist1 <= r1) & (dist2 > r2);

    case 'circular-annular-sector'

        % ...................... decode variables .........................

        center=domain_structure.center;
        r1=domain_structure.radii(1);
        r2=domain_structure.radii(2);
        alpha=domain_structure.angles(1);
        beta=domain_structure.angles(2);

        % ...................... indomain code ............................

        Xref=pts(:,1)-center(1); Yref=pts(:,2)-center(2);
        [TH,R] = cart2pol(Xref,Yref);

        inR=(R <= max(r1,r2)) & (R >= min(r1,r2));

        if not((abs(alpha) <= pi)  & (abs(beta) <= pi))
            error('alpha >= beta or at least one not in [-pi,pi]');
        end

        if (alpha < beta)
            % regular situation
            in_TH=(TH >= alpha) & (TH <= beta);
        else
            in_TH=not((TH >= beta) & (TH <= alpha));
        end

        in=inR & in_TH;



    case 'sector'

        % ...................... decode variables .........................

        center=domain_structure.center;
        r=domain_structure.radius;
        alpha=domain_structure.angles(1);
        beta=domain_structure.angles(2);

        % ...................... indomain code ............................

        Xref=pts(:,1)-center(1); Yref=pts(:,2)-center(2);
        [TH,R] = cart2pol(Xref,Yref);

        inR=(R <=r);

        if not((abs(alpha) <= pi)  & (abs(beta) <= pi))
            error('alpha >= beta or at least one not in [-pi,pi]');
        end

        if (alpha < beta)
            % regular situation
            in_TH=(TH >= alpha) & (TH <= beta);
        else
            in_TH=not((TH >= beta) & (TH <= alpha));
        end

        in=(R <=r) & in_TH;



    case  'asymmetric-annulus'

        % ...................... decode variables .........................

        centers=domain_structure.centers;
        rV=domain_structure.radii;

        % ...................... indomain code ............................

        C=centers(1,:);
        V=centers(2,:);
        r1=rV(1); r2=rV(2);

        X=pts(:,1); Y=pts(:,2);
        inCr1=(vecnorm([X-C(:,1) Y-C(:,2)],2,2) <= r1);
        inCr2=(vecnorm([X-V(:,1) Y-V(:,2)],2,2) <= r2);

        in=inCr1 & not(inCr2);

    case 'vertical-circular-zone'

        % ...................... decode variables .........................

        center=domain_structure.center;
        r=domain_structure.radius;
        alpha=domain_structure.angles(1);
        beta=domain_structure.angles(2);

        % ...................... indomain code ............................

        xmin=r*cos(beta); xmax=r*cos(alpha);

        Xref=pts(:,1)-center(1); Yref=pts(:,2)-center(2);
        [TH,R] = cart2pol(Xref,Yref);

        in=(R <= r) & ( Xref >= xmin ) & ( Xref <= xmax );

    case 'horizontal-circular-zone'

        % ...................... decode variables .........................

        center=domain_structure.center;
        r=domain_structure.radius;
        alpha=domain_structure.angles(1);
        beta=domain_structure.angles(2);

        % ...................... indomain code ............................

        ymin=r*sin(alpha); ymax=r*sin(beta);

        Xref=pts(:,1)-center(1); Yref=pts(:,2)-center(2);
        [TH,R] = cart2pol(Xref,Yref);

        in=(R <= r) & ( Yref >= ymin ) & ( Yref <= ymax );



    case 'symmetric-lens'

        % ...................... decode variables .........................

        a=domain_structure.center; a=abs(a);
        r=domain_structure.radius;

        center1=[-a 0];
        center2=[a 0];

        % ...................... indomain code ............................

        dist1=sqrt((pts(:,1)-center1(1)).^2+(pts(:,2)-center1(2)).^2);
        dist2=sqrt((pts(:,1)-center2(1)).^2+(pts(:,2)-center2(2)).^2);

        in=(dist1 <= r) & (dist2 <= r);


    case 'candy'

        % ...................... decode variables .........................

        a=domain_structure.center;
        r=domain_structure.radius;
        alpha=domain_structure.angle;

        % ...................... indomain code ............................

        error('Indomain routine for candy is not implemented');

    case {'NURBS','butterfly','asymmetric-circular-sector',...
            'circular-segment'}

        % ...................... indomain code ............................

        structure_RS=domain_structure.NURBS;
        [in,on]=indomain_NURBS(pts,structure_RS);

    case 'union-disks'

        % ...................... decode variables .........................

        centers=domain_structure.centers;
        rs=domain_structure.radii;

        % ...................... indomain code ............................

        Xc=centers(:,1); Yc=centers(:,2);
        X=pts(:,1); Y=pts(:,2);
        in=zeros(length(X),1);

        for k=1:length(rs)

            XLc=Xc(k); YLc=Yc(k); rL=rs(k);

            iout=find(in == 0);
            XL=X(iout); YL=Y(iout);
            r = vecnorm([XL-XLc YL-YLc],2,2);

            in(iout)=(r <= rL);

        end

        
    case 'polygcirc'

        % ...................... decode variables .........................

        ab=domain_structure.extrema;
        a=ab(1,:); b=ab(2,:);
        v=domain_structure.vertices;
        cc=domain_structure.center;
        r=domain_structure.radius;
        conv=domain_structure.convex;

        % ...................... indomain code ............................

        vertices=[a; v; b];

        X=pts(:,1); Y=pts(:,2);
        inP=inpolygon(X,Y,vertices(:,1),vertices(:,2));

        inD=( (X-cc(1)).^2+(Y-cc(2)).^2 <= r^2);

        if conv
            [THa,Ra] = cart2pol(a(1)-cc(1),a(2)-cc(2));
            [THb,Rb] = cart2pol(b(1)-cc(1),b(2)-cc(2));

            [TH,R] = cart2pol(X-cc(1),Y-cc(2));

            if THa <= THb
                inD=(R <= r) & not(TH >= THa & TH <= THb);
            else
                inD=(R <= r) & (TH >= THb & TH <= THa);
            end

            in= (inD | inP);

        else
            in= (not(inD) & inP);
        end



    otherwise
        error('indomain_routine not implemented for this domain');

end

