%% Setup
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

% Just 1 run for now..
% Variables
process_sig = 1e-3; % process noise variance
measure_sig = 1e-1; % measurement noise variance
Q = diag([0,0,0,process_sig^2, process_sig^2, process_sig^2]); % pn cov
R = diag([measure_sig^2, measure_sig^2, measure_sig^2]); % mn cov

% initial measurement (from initial truth)
h0 = range_angle_measurement(x0_tgt, x0_obs, t(1), R, mu);

% initial conditions
% Ukf cov
pcov = diag([1e-4, 1e-4, 1e-4, 1., 1., 1.]);
P = zeros(6,m); P(:,1) = diag(pcov); % cov storage

% ekf cov
pcov_ekf = pcov;
P_ekf = zeros(6,m); P_ekf(:,1)=diag(pcov); % ekf cov storage

% initial estimate
xhat0 = x0_tgt + 1e-3*randn(6,1);

% Unscented filter parameters
alpha=1e-3; beta=2; kappa=0; n=6; lam=alpha^2*(n+kappa)-n;

% Sigma weights
W0m = lam/(n+lam);
W0c = W0m + 1-alpha^2+beta;
Wim = 1/(2*n+2*lam);
Wic = Wim;
yez = zeros(3,2*n); % sigma point observation storage

% Generate truth and measurements
f = @(t, x) crtbp_natural_eom(t, x, mu, Q);
f_obs = @(t,x) crtbp_natural_eom(t, x, mu, diag(zeros(1,6))); 
X = zeros(n, m); X(:,1) = x0_tgt;
X_obs = zeros(n, m); X_obs(:,1) = x0_obs;
ym = zeros(3,m); ym(:,1) = h0;
for i=1:m-1
    % Truth (target)
    X(:,i+1) = rk4(f, X(:,i), t(i), dt);
    
    % Measurement
    X_obs(:,i+1) = rk4(f_obs, X_obs(:,i), t(i), dt); % propagate measurement
    ym(:,i+1) = range_angle_measurement(X(:,i+1), X_obs(:,i+1), t(i), R, mu);
end

% UF and EKF
Xhat = zeros(n, m); Xhat(:,1) = xhat0;
Xhat_ekf = zeros(n, m); Xhat_ekf(:,1) = xhat0;

% form augmented covariance

sigmat = zeros(12,25);
for i=1:m-1
    % UKF
    % Form augmented
    xa = [xhat0; zeros(6,1)];
    Pa = blkdiag(pcov, Q(4:6,4:6), R); % no process noise in Q(1:3) we hack this later

    % Cov decomp & sigma points
    gam = sqrt(n+lam);
    psq = gam*chol(Pa)';

    sig0 = xa;
    sigmat0 = real([xa+psq, xa-psq]);

    % propagate
    sigmat = [xa sigmat0];
    sigma_xs = sigmat(1:n,:); % xs only
    sigma_ws = sigmat(7:9,:); % ws opnly
    sigma_vs = sigmat(10:12,:);
    
    sigmatxs_ip1 = zeros(6,25); % only care about
    for ii=1:size(sigmat,2)
        sigma_w = sigma_ws(:,ii);

        xxi = rk4(f_obs, sigma_xs(:,ii), t(i), dt); % use no noise func.
        sigmatxs_ip1(:,ii) = xxi+[0;0;0;sigma_w]; %add noise
    end
    


    % EKF
    % propagate
    Xhat_ekf(:,i+1) = rk4(f, Xhat_ekf(:,i), t(i), dt);

    % estimate output (observation)
    ye_ekf = range_angle_measurement(Xhat_ekf(:,i+1), X_obs(:,i), t(i), R, mu);

    % covariance propagation
    F = crtbp_jacobian(Xhat_ekf(:,i+1), mu); % get jacobian
    phi = c2d(F, zeros(6,1), dt);
    pcov_ekf = phi*pcov_ekf*phi';

    % update
    h_ekf = range_angle_jacobian(Xhat_ekf(:,i+1), X_obs(:,i+1), mu, t(i));
    gain_ekf=pcov_ekf*h_ekf'*inv(h_ekf*pcov_ekf*h_ekf'+R);
    pcov_ekf=(eye(n)-gain_ekf*h_ekf)*pcov_ekf;
    P_ekf(:,i+1)=diag(pcov_ekf);
    Xhat_ekf(:,i+1)=Xhat_ekf(:,i+1)+(gain_ekf*(ym(:,i+1)-ye_ekf));
end

% Calculate errors ukf
err = Xhat-X;
% 3sig bounds
sig3 = P.^0.5 * 3;

% Calculate errors ekf
err_ekf = Xhat_ekf-X;

% 3sig bounds
sig3_ekf = P_ekf.^0.5 * 3;

figure
subplot(3,1,1)
hold on
l1=plot(t, err(1,:));
l2 = plot(t, sig3(1,:), '-r', t, -sig3(1,:), 'r'); l2=l2(1);
l3 = plot(t, err_ekf(1,:), '-k');
l4 = plot(t, sig3_ekf(1,:), '--k', t, -sig3_ekf(1,:), '--k'); l4 = l4(1);
legend([l1, l2, l3, l4], 'UKF err.', 'UKF 3sig', 'EKF err.', 'EKF 3sig', 'Location', 'northoutside', 'Orientation','horizontal')
title("Position Error")
hold off

subplot(3,1,2)
hold on
plot(t, err(2,:));
plot(t, sig3(2,:), '-r')
plot(t, -sig3(2,:), 'r');
plot(t, err_ekf(2,:), '-k');
plot(t, sig3_ekf(2,:), '--k', t, -sig3_ekf(2,:), '--k');
hold off

subplot(3,1,3)
hold on
plot(t, err(3,:))
plot(t, sig3(3,:), '-r', t, -sig3(3,:), 'r')
hold off

figure
subplot(3,1,1)
hold on
title("Velocity Error")
plot(t, err(4,:))
plot(t, sig3(4,:), '-r', t, -sig3(4,:), 'r')
hold off

subplot(3,1,2)
hold on
plot(t, err(5,:))
plot(t, sig3(5,:), '-r', t, -sig3(5,:), 'r')
hold off

subplot(3,1,3)
hold on
plot(t, err(6,:))
plot(t, sig3(6,:), '-r', t, -sig3(6,:), 'r')
hold off

figure
hold on
plot3(Xhat(1,:), Xhat(2,:), Xhat(3,:))
plot3(Xhat_ekf(1,:), Xhat_ekf(2,:), Xhat_ekf(3,:))
plot3(X(1,:), X(2,:), X(3,:))
hold off




% figure
% subplot(3,1,1)
% hold on
% plot(t, err(1,:))
% plot(t, sig3(1,:), '-r', t, -sig3(1,:), 'r')
% title("Position Error")
% hold off
% 
% subplot(3,1,2)
% hold on
% plot(t, err(2,:))
% plot(t, sig3(2,:), '-r', t, -sig3(2,:), 'r')
% hold off
% 
% subplot(3,1,3)
% hold on
% plot(t, err(3,:))
% plot(t, sig3(3,:), '-r', t, -sig3(3,:), 'r')
% hold off
% 
% figure
% subplot(3,1,1)
% hold on
% title("Velocity Error")
% plot(t, err(4,:))
% plot(t, sig3(4,:), '-r', t, -sig3(4,:), 'r')
% hold off
% 
% subplot(3,1,2)
% hold on
% plot(t, err(5,:))
% plot(t, sig3(5,:), '-r', t, -sig3(5,:), 'r')
% hold off
% 
% subplot(3,1,3)
% hold on
% plot(t, err(6,:))
% plot(t, sig3(6,:), '-r', t, -sig3(6,:), 'r')
% hold off
% 
% figure
% hold on
% plot3(Xhat_ekf(1,:), Xhat_ekf(2,:), Xhat_ekf(3,:))
% plot3(X(1,:), X(2,:), X(3,:))
% hold off
