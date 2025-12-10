clear; close all
syms x y z vx vy vz mu real

X = [x; y; z; vx; vy; vz];
t=0.0;
dx = crtbp_natural_eom(t, X, mu);
dx = simplify(dx);

F = jacobian(dx, X);

ccode(F, 'file', 'symbolic_crtbp_partials.txt')

function dx = crtbp_natural_eom(t, x, mu)
    %CRTBP_EOM Circular-restricted three-body problem ODE. Assumes canonical 
    % n.d. units, uncontrolled (natural) dynamics only.
    %   IN:
    %       - x: 6x1 synodic cartesian state vector [rx, ry, rz, vx, vy, vz]
    %       - mu: CRTBP mass parameter m2/(m1+m2) (0.01215 normally)
    %       - r_12: Earth-Moon distance (1 DU)
    %   OUT:
    %       - dx: 6x1 dx/dt
    
    % allocate
    
    % pull quantities
    r = x(1:3);
    vx = x(4); vy=x(5); vz=x(6);
   
    
    % dynamics
    gr = g(r,mu); % partial of pseudopotential wrt r
    
    dx123 = [vx; vy; vz];
    dx4 = gr(1) + 2*vy;
    dx5 = gr(2) - 2*vx;
    dx6 = gr(3);
    dx = [dx123; dx4; dx5; dx6];
end

function G = g(r, mu)
    rx = r(1); ry = r(2); rz = r(3);
    r1 = sqrt((rx+mu)^2 + ry^2 + rz^2); % sc-earth dist
    r2 = sqrt((rx+mu-1)^2 + ry^2 + rz^2); % sc-moon dist


    G1 = rx - (1-mu)*(rx+mu)/r1^3 - mu*(rx+mu-1)/r2^3;
    G2 = ry - (1-mu)*ry/r1^3 - mu*ry/r2^3;
    G3 = -(1-mu)*rz/r1^3 - mu*rz/r2^3;
    G = [G1; G2; G3];
end