% Testing natural crtbp eom
clear; close all
% CRTBP Parameters
G = 6.673e-20;
m_1 = 5.9742e24;         % kg
m_2 = 7.3483e22;         % kg
r_12 = 384400.0;          % km
mu = m_2 / (m_1 + m_2); % n.d.

%% Test 1: No occlusions
%time
tspan = [0, pi/2];
dt = 1e-4;

% observer
x0_obs = [1.1423846032; 0.0; 0.1597054213; 0.0; -0.2224918027; 0.0];

x0_tgt = [1.0186592988052636E+0	-2.1505075615805342E-27	-1.7967210088475610E-1	8.7422243833437732E-14	-9.5814062038783648E-2	1.3141536561558831E-12]';

wsig = 1e-3;
Q = diag([0,0,0,wsig^2, wsig^2, wsig^2]);

vsig = 1e-3;
R = diag([vsig^2, vsig^2, vsig^2]);

fun=@(t, s) crtbp_natural_eom(t, s, mu, Q);

T = 0:dt:tspan(2);
H = zeros(3,length(T));
X_obs = zeros(6,1); X_tgt = zeros(6,1);
X_obs(1:6,1) = x0_obs; X_tgt(1:6,1) = x0_tgt;

for i=1:length(T)-1
    t = T(i);
    
    % current states
    x_obs = X_obs(1:6,i);
    x_tgt = X_tgt(1:6,i);

    % measure
    H(1:3,i) = range_angle_measurement(x_tgt, x_obs, t, R, mu);
 
    % propagate
    X_obs(1:6,i+1) = rk4(fun, x_obs, t, dt);
    X_tgt(1:6,i+1) = rk4(fun, x_tgt, t, dt);

end

% last measurement
H(1:3,end) = range_angle_measurement(X_tgt(1:6,end), X_obs(1:6,end), t, R, mu);

figure
subplot(3,1,1)
plot(T, H(1,:))
ylabel("Range (DU)")

subplot(3,1,2)
plot(T, H(2,:).*180/pi)
ylabel("Az. deg")

subplot(3,1,3)
plot(T, H(3,:).*180/pi)
ylabel("El. deg")


figure
hold on
plot3(X_obs(1,:), X_obs(2,:), X_obs(3,:))
plot3(X_obs(1,1), X_obs(2,1), X_obs(3,1), '-kx')
plot3(X_tgt(1,:), X_tgt(2,:), X_tgt(3,:))
plot3(X_tgt(1,1), X_tgt(2,1), X_tgt(3,1), "-ro")
plot3(1.0-mu, 0.0, 0.0, 'kd')
xlabel("x"); ylabel("y"); zlabel("z")
view(45,45)
hold off