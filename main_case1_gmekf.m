%% Setup
clear; close all

% Don't have code exs from class to lean on, so doing this separately until
% I can code this satisfactorily

% Supress annoying warning
%#ok<*MINV> 

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
tf = pi;
t = 0:dt:tf;
m = length(t);
n=6;

process_sig = 1e-6; % process noise variance
measure_sig = 1e-1; % measurement noise variance
Q = diag([0,0,0,process_sig^2, process_sig^2, process_sig^2]); % pn cov
R = diag([measure_sig^2, measure_sig^2, measure_sig^2]); % mn cov

% initial measurement (from initial truth)
h0 = range_angle_measurement(x0_tgt, x0_obs, t(1), R, mu);

% initial conditions
pcov = diag([1e-6, 1e-6, 1e-6, 1e-3, 1e-3, 1e-3]);
P = zeros(6,m); P(:,1) = diag(pcov); % cov storage

% initial estimate
xhat0 = mvnrnd(x0_tgt, pcov, 1)'; % get realization from distribution

% Generate truth and measurements
f = @(t, x) crtbp_natural_eom(t, x, mu, Q);
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

% First, we need to approximate the distribution described by x0_tgt and
% pcov

% How many components (additional means) to use?
L = 5;
kk = randi(6);
[xmeans, ps, ws] = splitting_library(x0_tgt, pcov, L, kk);

covs = zeros(1,n,L);
for k=1:L
    covs(:,:,k)=ps(:,k);
end
gm = gmdistribution(xmeans', covs, ws); % vectors need to be horizontal for some reason
gms = cell(1,m-1); gms{1}=gm; % gmm storage in cell arr.
% note length is L-1 because we dont really care about gm at end of process

% now, pull the 'true' initial estimate from gm
xhat0 = random(gm)'; % transpose because gm stores horizontally

Xhat_gmekf = zeros(n,m); Xhat_gmekf(:,1)=xhat0;

for i=1:m-1
    
    ws = gm.ComponentProportion; % get weights of each component
    %L=L; % never changes in "baseline" GMEKF
    
    Xobsip1 = X_obs(:,i+1); % observer state
    ymip1 = ym(:,i+1); % measurement
    
    % propagate each mean & cov & estimate output
    ye_gmekf = zeros(3,L); % estimate output for each mean
    for k=1:L
        xm = xmeans(:,k);
        xmeans(:,k) = rk4(f, xm, t(i), dt);
        pk = diag(covs(1,:,k));
        
        Fm = crtbp_jacobian(xm, mu); % jacobian from prev. estimate
        phi = c2d(Fm, zeros(6,1), dt);
        Pkp = phi*pk*phi' + Q;
        covs(1,:,k) = diag(Pkp)';
    
        % est out
        ye_gmekf(:,k) = range_angle_measurement(xm, Xobsip1, t(i+1), R, mu);

    end


    % update
    betas = zeros(1,L); % term that will define the new weights
    for k=1:L
        
        % Some terms
        xm = xmeans(:,k);
        
        Hm = range_angle_jacobian(xm, Xobsip1, mu, t(i+1)); % linearized jacobian

        cov = diag(covs(1,:,k)); % propagated covariance of kth element

        Wk = Hm*cov*Hm'+R; % Term repeated a few times

        Kk = cov*Hm'*inv(Wk);  % gain
        
       
        % mean update
        xmeans(:,k)=xm+Kk*(-ye_gmekf(:,k));

        % Cov update
        covup = cov - Kk*Hm*cov;
        covs(1,:,k) = diag(covup); % save cov

        betak = det(2*pi*Wk)^(-1/2)*exp(-0.5 *(ymip1-ye_gmekf(:,k))'*inv(Wk)*(ymip1-ye_gmekf(:,k)));
        betas(:,k) = betak;

    end

    % now, set new weights
    sm = sum(betas.*ws);

    ws = ws.*betas./sm;


    gm = gmdistribution(xmeans', covs, eps+ws); % adding eps to let the gmdistribution object do the work for us
    
    gms{i} = gm;

    Xhat_gmekf(:,i+1) = random(gm);
    
end

all_ws = zeros(L,m-1);
for i=1:m-1
    ws = gms{i}.ComponentProportion;
    all_ws(:,i) = ws';
end

figure
plot(t(1:end-1), all_ws)


errs = X-Xhat_gmekf;

figure
subplot(3,1,1)
plot(t, errs(1,:))

subplot(3,1,2)
plot(t, errs(2,:))

subplot(3,1,3)
plot(t, errs(2,:))

