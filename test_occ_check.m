% Testing occlusion check
clear; close all
% CRTBP Parameters
G = 6.673e-20;
m_1 = 5.9742e24;         % kg
m_2 = 7.3483e22;         % kg
r_12 = 384400.0;          % km
mu = m_2 / (m_1 + m_2); % n.d.
DU = r_12;
TU = 1.0 / sqrt((G * (m_1 + m_2)) / DU^3);

% obs
x0_obs =  [1.1423846032; 0.0; 0.1597054213; 0.0; -0.2224918027; 0.0];
%x0_obs  = [0.823385182, 0.0, -0.022134068331103266, 0.0, 0.13408802571933368, 0.0]';
% DRO target state
x0 = [9.7651362015832144E-1	-5.7647806541303165E-28	-2.2746379925699658E-32	9.2638853182074272E-13	1.0468285446328065E+0	4.1733999115209969E-31	]';
x0 =[9.7746771543739996E-1	-6.0983443483835278E-28	-5.8351069377078023E-33	-1.6299713745531340E-12	1.0923700477597205E+0	2.7759195651165447E-30]';

sig = 0;
wvar = [0;0;0; sig^2; sig^2; sig^2]; % model error of each component
Q = diag(wvar); % covariance

dt = 60/TU;
tspan = [0, pi/2];

fun=@(t, s) crtbp_natural_eom(t, s, mu, Q);

T = 0:dt:tspan(2); m = length(T);
X = zeros(6,m);
X(1:6,1) = x0;
X_obs = zeros(6, m); X_obs(1:6,1) = x0_obs;

occ = -1*ones(1,m); occ(1) = is_occluded_by_moon(x0, x0_obs, mu);
H = zeros(3,m); H(:,1) = range_angle_measurement(x0, x0_obs, T(1), zeros(3,3), mu);

for i=1:length(T)-1
    % propagate
    t = T(i);
    x = X(1:6,i);
    X(1:6,i+1) = rk4(fun, x, t, dt);
    X_obs(1:6,i+1) = rk4(fun, X_obs(1:6,i), t, dt);

    % check occlusion and measure
    occ(i+1) = is_occluded_by_moon(X(1:6,i+1), X_obs(1:6,i+1), mu);
    H(:,i+1) = range_angle_measurement(X(1:6,i+1), X_obs(1:6,i+1), T(i+1), zeros(3,3), mu);
end

locc = logical(occ);
T_occ = T(locc);
hdims = 3;
H_occ = zeros(3, length(T_occ));
for i=1:hdims
    Hi = H(i,:);
    H_occ(i,:) = Hi(locc);
end


figure
subplot(3,1,1)
plot(T, H(1,:), T_occ, H_occ(1,:), '.r')
ylabel("Range (DU)")

subplot(3,1,2)
plot(T, H(2,:).*180/pi, T_occ, H_occ(2,:)*180/pi, '.r')
ylabel("Az. deg")

subplot(3,1,3)
plot(T, H(3,:).*180/pi, T_occ, H_occ(3,:)*180/pi, '.r')
ylabel("El. deg")

figure
hold on
plot3(X(1,:), X(2,:), X(3,:))
plot3(X(1,1), X(2,1), X(3,1), 'xk')
plot3(X_obs(1,:), X_obs(2,:), X_obs(3,:))
plot3(X_obs(1,1), X_obs(2,1), X_obs(3,1), 'dk')
hold off