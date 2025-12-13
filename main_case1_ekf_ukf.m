%% Setup
clear; close all

% num trials?
ntrials = 1;

% Pick Noise Level
noiselevel = "high"; % "low" or "high"
if noiselevel == "high"
    process_sig = 1e-3; % process noise variance
    measure_sig = 1e-1; % measurement noise variance
else
    process_sig = 0; % process noise variance
    measure_sig = 1e-4; % measurement noise variance
end


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



% Save data in cell arrays
sig3_all = cell(ntrials,1); % covariances sigma bounds
sig3_ekf_all = cell(ntrials,1);
Xes_all = cell(ntrials,1); %state errs
Xes_ekf_all = cell(ntrials, 1);

% the observer and truth states will still be calculated each time bc i
% dont want to rearrange things, so plots will just use the one made in the
% last trial

for trial=1:ntrials
    fprintf("\nTrial %i begin\n", trial)
    % Variables
    Q = diag([0,0,0,process_sig^2, process_sig^2, process_sig^2]); % pn cov
    R = diag([measure_sig^2, measure_sig^2, measure_sig^2]); % mn cov
    
    % initial measurement (from initial truth)
    h0 = range_angle_measurement(x0_tgt, x0_obs, t(1), R, mu);
    
    % initial conditions
    % Ukf cov
    pcov = diag([1e-6, 1e-6, 1e-6, 1e-3, 1e-3, 1e-3]);
    P = zeros(6,m); P(:,1) = diag(pcov); % cov storage
    
    % ekf cov
    pcov_ekf = pcov;
    P_ekf = zeros(6,m); P_ekf(:,1)=diag(pcov); % ekf cov storage
    
    % initial estimate
    xhat0 = mvnrnd(x0_tgt, pcov, 1)'; % get realization from distribution
    
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
    f_nonoise = @(t,x) crtbp_natural_eom(t, x, mu, zeros(6,6));
    X = zeros(n, m); X(:,1) = x0_tgt;
    X_obs = zeros(n, m); X_obs(:,1) = x0_obs;
    ym = zeros(3,m); ym(:,1) = h0;
    for i=1:m-1
        % Truth (target)
        X(:,i+1) = rk4(f, X(:,i), t(i), dt);
        
        % Measurement
        X_obs(:,i+1) = rk4(f_nonoise, X_obs(:,i), t(i), dt); % propagate observer
        ym(:,i+1) = range_angle_measurement(X(:,i+1), X_obs(:,i+1), t(i), R, mu);
    end
    
    % UF and EKF
    Xhat = zeros(n, m); Xhat(:,1) = xhat0;
    Xhat_ekf = zeros(n, m); Xhat_ekf(:,1) = xhat0;
    for i=1:m-1
        % UKF
        % Cov decomp & sigma points
        psquare = chol(pcov)';
        sigv=real([sqrt(n+lam)*psquare -sqrt(n+lam)*psquare]);
        xx0=Xhat(:,i);
        xx=sigv+kron(Xhat(:,i),ones(1,2*n));
    
        xxnext = rk4(f_nonoise, [xx0 xx], t(i), dt);
        xx0 = xxnext(:,1);
        xx = xxnext(:,2:2*n+1);
    
        Xhat(:,i+1) = W0m*xx0 + Wim*sum(xx,2);
    
        % cov
        pp0 = W0c*(xx0-Xhat(:,i+1))*(xx0-Xhat(:,i+1))';
        pmat = xx-kron(Xhat(:,i+1), ones(1,2*n));
        pcov = pp0+Wim*(pmat*pmat');
    
        for j=1:2*n
            yez(:,j)=range_angle_measurement(xx(:,j), X_obs(:,i+1), t(i+1), zeros(3,3), mu);
        end
        ye0 = range_angle_measurement(xx0, X_obs(:,i+1), t(i+1), zeros(3,3), mu);
        ye = W0m*ye0+Wim*sum(yez,2);
    
        % pyy
        pyy0 = W0c*((ye0-ye)*(ye0-ye)');
        pyymat = yez-ye;
        pyy = pyy0+Wim*(pyymat*pyymat');
    
        % pxy
        pxy0 = W0c*(xx0-Xhat(:,i+1))*(ye0-ye)';
        pxy = pxy0+Wim*pmat*pyymat';
    
        % Innovations
        pvv = pyy+R; % measurement error is linear!
    
        % Gain and update
        gain = real(pxy*inv(pvv));
        pcov = pcov-gain*pvv*gain';
        P(:,i+1) = diag(pcov);
        Xhat(:,i+1)=Xhat(:,i+1)+ gain*(ym(:,i+1)-ye);
    
    
        % EKF
        % propagate
        Xhat_ekf(:,i+1) = rk4(f_nonoise, Xhat_ekf(:,i), t(i), dt);
    
        % estimate output (observation)
        ye_ekf = range_angle_measurement(Xhat_ekf(:,i+1), X_obs(:,i+1), t(i+1), zeros(3,3), mu);
    
        % covariance propagation
        F = crtbp_jacobian(Xhat_ekf(:,i), mu); % get jacobian
        phi = c2d(F, zeros(6,1), dt);
        pcov_ekf = phi*pcov_ekf*phi'+Q*dt;
    
        % update
        h_ekf = range_angle_jacobian(Xhat_ekf(:,i+1), X_obs(:,i+1), mu, t(i+1));
        gain_ekf=pcov_ekf*h_ekf'*inv(h_ekf*pcov_ekf*h_ekf'+R); %#ok<*MINV>
        pcov_ekf=(eye(n)-gain_ekf*h_ekf)*pcov_ekf;
        P_ekf(:,i+1)=diag(pcov_ekf);
        Xhat_ekf(:,i+1)=Xhat_ekf(:,i+1)+(gain_ekf*(ym(:,i+1)-ye_ekf));
    end
    
    % save things from this trial
    sig3_all{trial} = P.^0.5*3;
    sig3_ekf_all{trial}=P_ekf.^0.5*3;
    Xes_all{trial}=Xhat-X;
    Xes_ekf_all{trial}=Xhat_ekf-X;


end

% Sum things from trials
sig3 = cat(3,sig3_all{:}); sig3=sum(sig3,3)/ntrials;
sig3_ekf = cat(3,sig3_ekf_all{:}); sig3_ekf=sum(sig3_ekf,3)/ntrials;
err = cat(3, Xes_all{:}); err=sum(err,3)/ntrials;
err_ekf = cat(3,Xes_ekf_all{:}); err_ekf=sum(err_ekf,3)/ntrials;

figure
subplot(3,1,1)
hold on
l1=plot(t, err(1,:));
l2 = plot(t, sig3(1,:), '-r', t, -sig3(1,:), 'r'); l2=l2(1);
l3 = plot(t, err_ekf(1,:), '-k');
l4 = plot(t, sig3_ekf(1,:), '--k', t, -sig3_ekf(1,:), '--k'); l4 = l4(1);
legend([l1, l2, l3, l4], 'UKF err.', 'UKF 3sig', 'EKF err.', 'EKF 3sig', 'Location', 'northoutside', 'Orientation','horizontal')
title("Position Error")
ylabel("x (DU)")
hold off

subplot(3,1,2)
hold on
plot(t, err(2,:));
plot(t, sig3(2,:), '-r')
plot(t, -sig3(2,:), 'r');
plot(t, err_ekf(2,:), '-k');
plot(t, sig3_ekf(2,:), '--k', t, -sig3_ekf(2,:), '--k');
ylabel("y (DU)")
hold off

subplot(3,1,3)
hold on
plot(t, err(3,:))
plot(t, sig3(3,:), '-r', t, -sig3(3,:), 'r')
plot(t, err_ekf(3,:), '-k');
plot(t, sig3_ekf(3,:), '--k', t, -sig3_ekf(3,:), '--k');
ylabel("z (DU)")
hold off

figure
subplot(3,1,1)
hold on
title("Velocity Error")
l1 = plot(t, err(4,:)); 
l2 = plot(t, sig3(4,:), '-r', t, -sig3(4,:), 'r'); l2=l2(1);
l3 = plot(t, err_ekf(4,:), '-k');
l4 = plot(t, sig3_ekf(4,:), '--k', t, -sig3_ekf(4,:), '--k'); l4 = l4(1);
ylabel("vx (DU/TU)")
legend([l1, l2, l3, l4], 'UKF err.', 'UKF 3sig', 'EKF err.', 'EKF 3sig', 'Location', 'northoutside', 'Orientation','horizontal')
hold off

subplot(3,1,2)
hold on
plot(t, err(5,:))
plot(t, sig3(5,:), '-r', t, -sig3(5,:), 'r')
plot(t, err_ekf(5,:), '-k');
plot(t, sig3_ekf(5,:), '--k', t, -sig3_ekf(5,:), '--k');
ylabel("vy (DU/TU)")
hold off

subplot(3,1,3)
hold on
plot(t, err(6,:))
plot(t, sig3(6,:), '-r', t, -sig3(6,:), 'r')
plot(t, err_ekf(6,:), '-k');
plot(t, sig3_ekf(6,:), '--k', t, -sig3_ekf(6,:), '--k');
ylabel("vz (DU/TU)")
xlabel("Time (TU)")
hold off

figure
hold on
plot3(Xhat(1,:), Xhat(2,:), Xhat(3,:))
plot3(Xhat_ekf(1,:), Xhat_ekf(2,:), Xhat_ekf(3,:))
plot3(X(1,:), X(2,:), X(3,:))
legend('UKF', 'EKF', 'Truth')
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
