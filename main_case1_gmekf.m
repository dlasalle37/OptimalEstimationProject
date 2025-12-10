%% Setup
clear; close all

% Don't have code exs from class to lean on, so doing this separately until
% I can code this satisfactorily


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

process_sig = 1e-3; % process noise variance
measure_sig = 1e-1; % measurement noise variance
Q = diag([0,0,0,process_sig^2, process_sig^2, process_sig^2]); % pn cov
R = diag([measure_sig^2, measure_sig^2, measure_sig^2]); % mn cov

% initial measurement (from initial truth)
h0 = range_angle_measurement(x0_tgt, x0_obs, t(1), R, mu);

% initial conditions
pcov = diag([1e-4, 1e-4, 1e-4, 1., 1., 1.]);
P = zeros(6,m); P(:,1) = diag(pcov); % cov storage

% initial estimate
xhat0 = mvnrnd(x0_tgt, pcov, 1)'; % get realization from distribution

% Generate truth and measurements
f = @(t, x) crtbp_natural_eom(t, x, mu, Q);
f_nonoise = @(t,x) crtbp_natural_eom(t, x, mu, zeros(6,6));
f_obs = @(t,x) crtbp_natural_eom(t, x, mu, zeros(6,6)); 
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

% Initial GMEKF parameters