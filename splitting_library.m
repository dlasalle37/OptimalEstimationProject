function [ms, Ps, ws] = splitting_library(m, P, L, k)
%SPLIT_GAUSSIAN Split a distribution with mean m and cov P into L gaussians
% note that covariance of each split gaussian are equal
% based on "splitting libraries" as described in 
% Entropy-Based Approach for Uncertainty Propagation of Nonlinear Dynamical Systems
% by DeMars, Bishop, and Jah

% IN: 
%   m: mean (nx1)
%   P: covariance (nxn)
%   L: which library (3 or 5)
%   k: which eigenvalue/vector to use to split along
    
    % parameter values from splitting library in above work
    if L == 3 % beta = 2, lambda = 0.001
        alpha_tilde = [0.2252246249; 0.5495507502; 0.2252246249];
        m_tilde = [-1.0575154615; 0.0; 1.0575154615];
        sig_tilde = [ 0.6715662887; 0.6715662887; 0.6715662887];
    elseif L == 5 % beta=2, lambda =0.0025
        alpha_tilde = [0.0763216491; 0.2474417860; 0.3524731300; 0.2474417860; 0.0763216491];
        m_tilde = [-1.6899729111; -0.8009283834; 0; -0.8009283834; -1.6899729111];
        sig_tilde = [0.4422555386; 0.4422555386; 0.4422555386; 0.4422555386; 0.4422555386];
    else
        error("Only values of 3 or 5 are allowed for L!")
    end

    alpha = 1; % assuming weight of original component is 1 (this shouldnt cause loss of generality, just need to scale appropriately later)

    n = length(m); % dimension of mean vector
    
    ms = zeros(n, L);
    Ps = zeros(n, L); % diagonal variances only;
    ws = zeros(1,L); % weights

    % doing the sqrt method with eigen decomp as done in cited paper
    % Assuming P is symmetric (which should always be the case since its a
    % covariance):
    [V, LAM] = eig(P);
    lams = diag(LAM);
    lamk = lams(k);
    vk = V(:,k);
    

    for i = 1:L
        ws(i) = alpha_tilde(i)*alpha;
        ms(:,i) = m + sqrt(lamk)*m_tilde(i)*vk;
           
        % construct Lambda_i matrix
        for j = 1:n
            if j==k
                lams(j)=lams(j)*sig_tilde(i)^2;
            end
        end
        LAMi = diag(lams);

        Pi = V*LAMi*V';
        Ps(:,i) = diag(Pi);

    end

    mss{k} = ms;
    Pss{k} = Ps;
    wss{k}=ws;
end

