function xout = rk4(fun, x, t, dt)
%RK4 Summary Runge-kutta 4th order step
%   IN:
%       - fun(t, x): function to step

    k1 = fun(t, x);
    k2 = fun(t+0.5*dt, x+k1*dt/2);
    k3 = fun(t+0.5*dt, x+k2*dt/2);
    k4 = fun(t+dt, x+k3*dt);

    xout = x + (dt/6)*(k1+2*k2+2*k3+k4);


end

