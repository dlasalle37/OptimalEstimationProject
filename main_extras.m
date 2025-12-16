% Some extra figures
clear; close all
% CRTBP params
G = 6.673e-20;
m_1 = 5.9742e24;         % kg
m_2 = 7.3483e22;         % kg
r_12 = 384400.0;          % km
mu = m_2 / (m_1 + m_2); % n.d.
DU = r_12;
TU = 1.0 / sqrt((G * (m_1 + m_2)) / DU^3);
% SC Initial states
% Observer (L2 halo)
x0_obs = [1.1423846032; 0.0; 0.1597054213; 0.0; -0.2224918027; 0.0];

% Target (NRHO)
x0_tgt = [1.0186592988052636E+0	-2.1505075615805342E-27	-1.7967210088475610E-1	8.7422243833437732E-14	-9.5814062038783648E-2	1.3141536561558831E-12]';



% sampling and interval time
dt = 60/TU; % every minute
tf = pi/2;
t = 0:dt:tf;
m = length(t);
n=6;
R = zeros(3,3);
Q=zeros(6,6);
% initial measurement (from initial truth)
h0 = range_angle_measurement(x0_tgt, x0_obs, t(1), R, mu);
% Generate truth and measurements
f = @(t, x) crtbp_natural_eom(t, x, mu, Q);
f_nonoise = @(t,x) crtbp_natural_eom(t, x, mu, zeros(6,6)); 
X = zeros(n, m); X(:,1) = x0_tgt;
X_obs = zeros(n, m); X_obs(:,1) = x0_obs;
ym = zeros(3,m); ym(:,1) = h0;

for i=1:m-1
    % Truth (target)
    X(:,i+1) = rk4(f, X(:,i), t(i), dt);
    
    % Measurement
    X_obs(:,i+1) = rk4(f_nonoise, X_obs(:,i), t(i), dt); % propagate measurement
    ym(:,i+1) = range_angle_measurement(X(:,i+1), X_obs(:,i+1), t(i), R, mu);
end

figure
hold on
l1 = plot3(X(1,:), X(2,:), X(3,:));
plot3(x0_tgt(1), x0_tgt(2), x0_tgt(3), 'xk')
l2 = plot3(X_obs(1,:), X_obs(2,:), X_obs(3,:));
l3 = plot3(x0_obs(1), x0_obs(2), x0_obs(3), 'xk');
[Xx,Yy,Zz] = sphere(50); Xx = Xx*1738/DU+1.0-mu;
Yy = 1738/DU*Yy; Zz = 1738/DU*Zz;
surf(Xx,Yy,Zz); axis equal
l4 = plot3(1.0-mu, 0.0, 0.0, '.k');
legend([l1, l2, l3, l4], ["Target", "Observer", "Initial Pos.", "Moon"])
xlabel("X (DU)"); ylabel("Y (DU)"); zlabel("Z (DU)")
grid on
hold off


% Observer (L2 Halo)
x0_obs =  [1.108176760562800, -0.117270070982379, 0.106310036921042, -0.096157733123961, -0.115664544723324, -0.166026773608929]';

% Target (dro) % jc = 3.996
x0_tgt = [9.7651362015832144E-1	-5.7647806541303165E-28	-2.2746379925699658E-32	9.2638853182074272E-13	1.0468285446328065E+0	4.1733999115209969E-31	]';

% sampling and interval time
dt = 60/TU; % every minute
tf = pi/4;
t = 0:dt:tf;
m = length(t);
n=6;
X = zeros(n, m); X(:,1) = x0_tgt;
X_obs = zeros(n, m); X_obs(:,1) = x0_obs;
ym = zeros(3,m); ym(:,1) = h0;
for i=1:m-1
    % Truth (target)
    X(:,i+1) = rk4(f, X(:,i), t(i), dt);
    
    % Measurement
    X_obs(:,i+1) = rk4(f_nonoise, X_obs(:,i), t(i), dt); % propagate measurement
    ym(:,i+1) = range_angle_measurement(X(:,i+1), X_obs(:,i+1), t(i), R, mu);
end

figure
hold on
l1 = plot3(X(1,:), X(2,:), X(3,:));
plot3(x0_tgt(1), x0_tgt(2), x0_tgt(3), 'xk')
l2 = plot3(X_obs(1,:), X_obs(2,:), X_obs(3,:));
l3 = plot3(x0_obs(1), x0_obs(2), x0_obs(3), 'xk');
[Xx,Yy,Zz] = sphere(50); Xx = Xx*1738/DU+1.0-mu;
Yy = 1738/DU*Yy; Zz = 1738/DU*Zz;
surf(Xx,Yy,Zz); axis equal
l4 = plot3(1.0-mu, 0.0, 0.0, '.k');
legend([l1, l2, l3, l4], ["Target", "Observer", "Initial Pos.", "Moon"])
xlabel("X (DU)"); ylabel("Y (DU)"); zlabel("Z (DU)")
grid on
hold off
