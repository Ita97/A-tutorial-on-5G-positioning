function [x,y] = getRSTDCurve(gNB1,gNB2,rstd)
%   [X,Y] = getRSTDCurve(GNB1,GNB2,RSTD) returns the x- and y-coordinates
%   of a hyperbola equation corresponding to an RSTD value, given the pair
%   of gNB positions gNB1 and gNB2.

    % Calculate vector between two gNBs
    delta = gNB1 - gNB2;

    % Express distance vector in polar form
    [phi,r] = cart2pol(delta(1),delta(2));
    rd = (r+rstd)/2;

    % Compute the hyperbola parameters
    a = (r/2)-rd;         % Vertex
    c = r/2;              % Focus
    b = sqrt(c^2-a^2);    % Co-vertex
    hk = (gNB1 + gNB2)/2;
    mu = -2:1e-3:2;

    % Get x- and y- coordinates of hyperbola equation
    x = (a*cosh(mu)*cos(phi)-b*sinh(mu)*sin(phi)) + hk(1);
    y = (a*cosh(mu)*sin(phi)+b*sinh(mu)*cos(phi)) + hk(2);

end