function [xhat, Phat] = produce_state_estimate(xms, ps, ws)
% From DeMars GMEKF paper    

    [n, L] = size(xms);
    
    % means
    xhat = zeros(n,1);
    Phat = zeros(n,n);
    for k=1:L
        wm = ws(k);
        xm = xms(:,k);
        xhat = xhat + wm*xm;
        
        Pm = diag(ps(:,k));
        Phat = Phat + wm*Pm + wm*xm*xm';
    end
    Phat = Phat - xhat*xhat';

    
end