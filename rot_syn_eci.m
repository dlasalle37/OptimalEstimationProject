function x_eci = rot_syn_eci(x_syn, t, mu)
%ROT_SYN_ECI Rotate state `x` to eci frame at time `t`. Assumes initially
%aligned at t=0.

% since im using canonical units, this is easy. TU is defined such that a
% full rotation completes every 2*pi TU
% first, find rotation angle
rotation_angle = mod(t, 2*pi);
rotation_angle=t;

r_syn = x_syn(1:3);
v_syn = x_syn(4:6);

sw = sin(rotation_angle);
cw = cos(rotation_angle); 
 
r_eci = [cw -sw 0; sw cw 0; 0 0 1]*(r_syn+[mu;0;0]); % rotation matrix and accounting for shift 
v_eci = v_syn + cross([0;0;1], r_syn); % transport theorem

x_eci = [r_eci; v_eci];

end

