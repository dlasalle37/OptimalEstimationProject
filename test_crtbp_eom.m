% Testing natural crtbp eom
clear; close all
% CRTBP Parameters
G = 6.673e-20;
m_1 = 5.9742e24;         % kg
m_2 = 7.3483e22;         % kg
r_12 = 384400.0;          % km
mu = m_2 / (m_1 + m_2); % n.d.
DU = r_12;
TU = 1.0 / sqrt((G * (m_1 + m_2)) / DU^3);

% NRHO State (parrish)
r0 = [0.9956461791199591; -0.04622742816025321; -0.05094004418576085];
v0 = [-0.08748056716039979; 0.11304919197855198; 0.4906469979990478];
x0 = [r0;v0];

% another nrho state (6.50 period, should be approx 9:2, JC)
x0 = [1.0186592988052636E+0	-2.1505075615805342E-27	-1.7967210088475610E-1	8.7422243833437732E-14	-9.5814062038783648E-2	1.3141536561558831E-12]';

% L1 Halo state
x0 = [0.823385182, 0.0, -0.022134068331103266, 0.0, 0.13408802571933368, 0.0];

% DRO State
 x0 = [1.0773094648; 0.0; 0.0; 0.0; -0.4697376289; 0.0];

sig = 0;
wvar = [0;0;0; sig^2; sig^2; sig^2]; % model error of each component
Q = diag(wvar); % covariance

dt = 60/TU;
tspan = [0, pi/2];

fun=@(t, s) crtbp_natural_eom(t, s, mu, Q);

T = 0:dt:tspan(2);
X = zeros(6,1);
X(1:6,1) = x0;
for i=1:length(T)-1
    t = T(i);
    x = X(1:6,i);
    X(1:6,i+1) = rk4(fun, x, t, dt);
end

[T_45, X_45] = ode45(fun, [tspan(1):dt:tspan(2)], x0, odeset('RelTol',1e-12, 'AbsTol',1e-12));
X_45 = X_45';
T_45=T_45';

figure
hold on
plot3(X(1,:), X(2,:), X(3,:))
plot3(X_45(1,:), X_45(2,:), X_45(3,:), '--r')
plot3(1.0-mu, 0.0, 0.0, 'xk')
xlabel('X'); ylabel('Y'); zlabel("Z")
view(34,34)

figure

subplot(3,1,1)
plot(T, X(1,:)-X_45(1,:))
title_str = sprintf("ode45 vs rk4, dt=%4.1s", dt);
title(title_str)
subplot(3,1,2)
plot(T, X(2,:)-X_45(2,:))
subplot(3,1,3)
plot(T, X(3,:)-X_45(3,:))

% test the syn->eci conversion
rx_eci = zeros(1,length(T_45)); ry_eci = zeros(1,length(T_45)); rz_eci = zeros(1,length(T_45));
for i = 1:length(T_45)
    t=T_45(i);
    x_syn = X_45(:,i);
    x_eci = rot_syn_eci(x_syn, t, mu);
    rx_eci(i) = x_eci(1);
    ry_eci(i) = x_eci(2);
    rz_eci(i) = x_eci(3);
end

figure
plot3(rx_eci, ry_eci, rz_eci)