function dx = crtbp_natural_eom(t, x, mu, Q)
    %CRTBP_EOM Circular-restricted three-body problem ODE. Assumes canonical 
    % n.d. units, uncontrolled (natural) dynamics only.
    %   IN:
    %       - x: 6x1 synodic cartesian state vector [rx, ry, rz, vx, vy, vz]
    %       - mu: CRTBP mass parameter m2/(m1+m2) (0.01215 normally)
    %       - r_12: Earth-Moon distance (1 DU)
    %       - Q: 6x6 zero-mean gaussian noise covariance
    %   OUT:
    %       - dx: 6x1 dx/dt
    
    % allocate
    [m, n] = size(x);
    dx = zeros(m,n);
    
    % loop through columns of x
    for i = 1:n
        % pull quantities
        r = x(1:3,i);
        vx = x(4,i); vy=x(5,i); vz=x(6,i);
        
        % Form noise vec
        w = Q*randn(m,1);
        
        % dynamics
        gr = g(r,mu); % partial of pseudopotential wrt r
        
        dx(1:3,i) = [vx; vy; vz];
        dx(4,i) = gr(1) + 2*vy;
        dx(5,i) = gr(2) - 2*vx;
        dx(6,i) = gr(3);
    
        % add the noise
        dx(:,i) = dx(:,i)+w;
    end
end

function G = g(r, mu)
    G = zeros(3,1);
    rx = r(1); ry = r(2); rz = r(3);
    r1 = sqrt((rx+mu)^2 + ry^2 + rz^2); % sc-earth dist
    r2 = sqrt((rx+mu-1)^2 + ry^2 + rz^2); % sc-moon dist


    G(1) = rx - (1-mu)*(rx+mu)/r1^3 - mu*(rx+mu-1)/r2^3;
    G(2) = ry - (1-mu)*ry/r1^3 - mu*ry/r2^3;
    G(3) = -(1-mu)*rz/r1^3 - mu*rz/r2^3;
end