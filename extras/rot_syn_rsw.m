function R_ECI_RSW = rot_syn_rsw(x_syn,t, mu)
%ROT_SYN_RSW rotate state `x` to rsw frame at time `t`. Return rotation
%matrix
%  Assumes ECI frame is aligned with synodic at t=0.
    
    % perform first rotation
    x_eci = rot_syn_eci(x_syn, t, mu);
    
    % partition
    r_eci = x_eci(1:3);
    v_eci = x_eci(4:6);
    
    % R basis
    R = r_eci/norm(r_eci);
    
    % W basis
    cprod = cross(r_eci, v_eci);
    W = cprod/norm(cprod);

    % S Basis
    S = cross(W, R);
    
    % Construct Matrix
    R_ECI_RSW = [R S W];

end

