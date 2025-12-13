function [xhat, Phat] = produce_state_estimate(xms, ps, ws)
% From DeMars GMEKF paper  
% produce a state estimate from a GMM, or single gaussian 
% IN: 
%   xms: nxL matrix of means. n is state size, L is number of mixture
%   elements
%   

    [n, L] = size(xms);
    
    % means
    xhat = zeros(n,1);
    Phat = zeros(n,n);
    for k=1:L
        wm = ws(k);
        xm = xms(:,k);
        xhat = xhat + wm*xm;
        
        %Pm = diag(ps(:,k));
        Pm = ps{k};
        Phat = Phat + wm*Pm + wm*(xm*xm');
    end
    Phat = Phat - xhat*xhat';

    
end