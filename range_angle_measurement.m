function h = range_angle_measurement(x, x_obs, t, R, mu)
%RANGE_ANGLE_MEASUREMENT return 1x3 relative range and angle measurements
%from obs. This is the measurement model
%   IN:
%       x: target state
%       x_obs: observer state
%       t: time
%       R: 3x3 zero-mean gaussian measurement noise covariance

r = x(1:3);
r_obs = x_obs(1:3);
    
% get rotation matrix
R_ECI_RSW = rot_syn_rsw(x_obs, t, mu);

% get range vector
rho_syn = r-r_obs;
rho_rsw = R_ECI_RSW*rho_syn;
rho_r = rho_rsw(1); rho_s=rho_rsw(2); rho_w=rho_rsw(3);

% get angles and range
rho_mag = norm(rho_rsw);
h = [rho_mag; atan(rho_s/rho_r); asin(rho_w/rho_mag)]; %[range, az., el.]'

% add noise
v = R*randn(3,1);
h = h + v;

end

