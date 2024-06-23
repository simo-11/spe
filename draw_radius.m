function [n,points]=draw_radius(pt,r,theta,n,ccw,phi)
    % based on draw_radius in utils.py of section-properties
    arguments
        pt % start point
        r % radius
        theta % start angle in radians
        n % number of points to generate
        ccw=true % counter clockwise i.e. angle increases
        phi=0.5*pi % angle to cover in radians
    end
    if r==0
        n=1;
        points=pt;
        return;
    end
    if ccw
        mult = 1;
    else
        mult = -1;
    end
    points=zeros(n,2);
    % calculate radius of points
    for i=1:n
        % determine angle
        t = theta + mult * (i-1) * 1.0 / max(1, n - 1) * phi;
        points(i,1) = pt(1) + r * cos(t);
        points(i,2) = pt(2) + r * sin(t);
    end
end