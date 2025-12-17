%% Setup
clear; close all

% This filter is a little more involved, so coding it separately

% Supress annoying warning
%#ok<*MINV> 

% num trials?
ntrials = 100;

% setting to true will overwrite plots/
save_figs = false;


% Pick Noise Level
noiselevel = "high"; % "low" or "high"
if noiselevel == "high"
    process_sig = 1e-4; % sqrt process noise variance
    measure_sig = 1e-2; % sqrt measurement noise variance
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
n=6;

% Saving
Xes_all = cell(1,ntrials); % estimate errors for all
sig3_all = cell(1,ntrials); % sig3 for all
gmms_all = cell(1,ntrials); % mixture models for all

for jj=1:ntrials
    fprintf("Beginning Trial %i\n", jj)
    Q = diag([0,0,0,process_sig^2, process_sig^2, process_sig^2]); % pn cov
    R = diag([measure_sig^2, measure_sig^2, measure_sig^2]); % mn cov
    
    % initial measurement (from initial truth)
    h0 = range_angle_measurement(x0_tgt, x0_obs, t(1), R, mu);
    
    % initial conditions
    pcov = diag([1e-6, 1e-6, 1e-6, 1e-3, 1e-3, 1e-3]);
    P = zeros(6,m); P(:,1) = diag(pcov); % cov storage
    
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
    
    % Initial GMEKF parameters
    
    % First, we need to approximate the distribution described by x0_tgt and
    % pcov
    
    % How many components (additional means) to use?
    xhat0 = x0_tgt+chol(pcov)'*randn(6,1);
    L = 5;
    kk = randi(6);
    [xmeans, ps, ws] = splitting_library(xhat0, pcov, L, kk);
    
    ps2 = cell(1,L);
    for i=1:L
        ps2{i}=diag(ps(:,i));
    end
    ps = ps2;
    
    % Here's how I'm going to organize the GMMs
    % the for a given trial, the GMMs will be a cell(3,m-1);
    % each column will consist of 3 matrices in each cell
    % cell in row 1: nxL matrix of means
    % cell in row 2: 1xL cell arr. of covariances
    % cell in row 3: 1xL matrix of weights
    
    gms = cell(3,m-1);
    gms{1,1} = xmeans;
    gms{2,1} = ps;
    gms{3,1} = ws;
    
    % now, pull the 'true' initial estimate from gmm
    [xhat0, Phat0] = produce_state_estimate(gms{1,1}, gms{2,1}, gms{3,1}); 
    Xhat_gmekf = zeros(n,m); Xhat_gmekf(:,1)=xhat0;
    P_combined_gmekf = zeros(n,m); P_combined_gmekf(:,1)=diag(Phat0);
    for i=1:m-1
        
        ws = gms{3,i};
        
        Xobsip1 = X_obs(:,i+1); % observer state
        ymip1 = ym(:,i+1); % measurement
        
        % propagate each mean & cov & estimate output
        ye_gmekf = zeros(3,L); % estimate output for each mean
        for k=1:L
            xm = xmeans(:,k); % xm_(i)^+ (previous step posterior)
            xmeans(:,k) = rk4(f_nonoise, xm, t(i), dt); % xm(i+1)^- (current step priori)
            %pk = diag(ps(:,k)); % previous step posterior
            pk = ps{k};
            Fm = crtbp_jacobian(xm, mu); % jacobian from prev. mean
            phi = c2d(Fm, zeros(6,1), dt);
            Pkp = phi*pk*phi' + Q*dt;
            %ps(:,k) = diag(Pkp); % current step priori
            ps{k}=Pkp;
            % est out
            ye_gmekf(:,k) = range_angle_measurement(xmeans(:,k), Xobsip1, t(i+1), zeros(3,3), mu);
    
        end
    
    
        % update
        betas = zeros(1,L); % term that will define the new weights
        for k=1:L
            
            % Some terms
            xm = xmeans(:,k); % prediction mean (mi^-)
            Hm = range_angle_jacobian(xm, Xobsip1, mu, t(i+1)); % jacobian
            %cov = diag(ps(:,k)); % propagated covariance of kth element
            cov = ps{k};
            Wk = Hm*cov*Hm'+R; % Term repeated a few times
            
            % gain
            Kk = cov*Hm'*inv(Wk);    
           
            % mean update
            xmeans(:,k)=xm+Kk*(ymip1-ye_gmekf(:,k));
    
            % Cov update
            cov = (eye(n)-Kk*Hm)*cov;
            %ps(:,k) = diag(cov); % save cov
            ps{k}=cov;
            betak = det(2*pi*Wk)^(-1/2)*exp(-0.5 *(ymip1-ye_gmekf(:,k))'*inv(Wk)*(ymip1-ye_gmekf(:,k)));
            betas(:,k) = betak;
    
        end
    
        % now, set new weights
        sm = sum(betas.*ws);
    
        ws = ws.*betas./sm;
        
        gms{1,i+1} = xmeans;
        gms{2,i+1} = ps;
        gms{3,i+1} = ws;
    
        [xh, ph] = produce_state_estimate(xmeans, ps, ws);
        Xhat_gmekf(:,i+1) = xh;
        P_combined_gmekf(:,i+1) = diag(ph);
    end

    % save trial data
    Xes_all{jj} = Xhat_gmekf-X;
    sig3_all{jj} = P_combined_gmekf.^0.5*3;
    gmms_all{jj} = gms;

end

% Plot example weights from last trial
all_ws = zeros(L,m-1);
for i=1:m-1
    ws = gms{3,i};
    all_ws(:,i) = ws';
end

% average weights
ws_avg = zeros(5,m);
for jj=1:ntrials
    gm = gmms_all{jj};
    for i = 1:m
        ws_avg(:,i) = ws_avg(:,i)+gm{3,i}';
    end
end
ws_avg = ws_avg/ntrials;
fwavg = figure;
plot(t,ws_avg)
xlabel("Time (TU)")
legend("w_1", "w_2", "w_3", "w_4", "w_5", 'Location','eastoutside')
ylabel("w_i")

fwex = figure;
plot(t(1:end-1), all_ws)
xlabel("Time (TU)")
legend("w_1", "w_2", "w_3", "w_4", "w_5", 'Location','eastoutside')
ylabel("w_i")

% grab average errors and 3sig bounds
sig3 = cat(3,sig3_all{:}); sig3=sum(sig3,3)/ntrials;
errs = cat(3, Xes_all{:}); errs=sum(errs,3)/ntrials;

f1 = figure;
subplot(3,1,1)
hold on
l1=plot(t, errs(1,:));
l2=plot(t, sig3(1,:), '-r');
plot(t, -sig3(1,:), '-r')
legend([l1,l2], ["Error", "Sig3"], 'Orientation','horizontal', 'Location','northoutside')
ylabel("X (DU)")
hold off

subplot(3,1,2)
hold on
plot(t, errs(2,:))
plot(t, sig3(2,:), '-r')
plot(t, -sig3(2,:), '-r')
ylabel("Y (DU)")
hold off

subplot(3,1,3)
hold on
plot(t, errs(3,:))
plot(t, sig3(3,:), '-r')
plot(t, -sig3(3,:), '-r')
ylabel("Z (DU)")
xlabel("Time (TU)")
hold off

f2 = figure;

subplot(3,1,1)
hold on
l1=plot(t, errs(4,:));
l2=plot(t, sig3(4,:), '-r');
plot(t, -sig3(4,:), '-r')
ylabel("X (DU/TU)")
legend([l1,l2], ["Error", "Sig3"], 'Orientation','horizontal', 'Location','northoutside')
hold off

subplot(3,1,2)
hold on
plot(t, errs(5,:))
plot(t, sig3(5,:), '-r')
plot(t, -sig3(5,:), '-r')
ylabel("Y (DU/TU)")
hold off

subplot(3,1,3)
hold on 
plot(t, errs(6,:))
plot(t, sig3(6,:), '-r')
plot(t, -sig3(6,:), '-r')
ylabel("Z (DU/TU)")
xlabel("Time (TU)")
hold off


f3 = figure;
hold on
plot3(Xhat_gmekf(1,:), Xhat_gmekf(2,:), Xhat_gmekf(3,:))
plot3(X(1,:), X(2,:), X(3,:))
legend('GMEKF', 'Truth')
hold off

f4 = figure;
subplot(3,1,1)
plot(t, ym(1,:))
ylabel("Range (DU)")

subplot(3,1,2)
plot(t,ym(2,:)*180/pi)
ylabel("Az. (deg)")

subplot(3,1,3)
plot(t,ym(3,:)*180/pi)
ylabel("El. (deg)")
xlabel("Time (TU)")

if save_figs == true
exportgraphics(f1, "plots/case_1_"+noiselevel+"_mc_gmekf.png", 'Resolution',600)
    exportgraphics(f2, "plots/case_1_"+noiselevel+"_mc_gmekf_vel.png", 'Resolution', 600)
    exportgraphics(f3, "plots/case_1_"+noiselevel+"_traj_gmekf.png", 'Resolution', 600)
    exportgraphics(f4, "plots/case_1_measurement_"+noiselevel+".png", 'Resolution', 600)
    exportgraphics(fwex, "plots/case_1_"+noiselevel+"_ws.png", 'Resolution', 600)
end