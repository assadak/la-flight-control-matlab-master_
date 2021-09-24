function [x_next] = RK4(x,u,p, h,f)
%
% Inputs : 
%    x, u, p current state and input
%    h    sample period
%    f    continuous time dynamics f(x,u,p)
% Returns
%    State h seconds in the future
%

% Runge-Kutta 4 integration
% write your function here
   k1 = f(x,        u, p);
   k2 = f(x+h/2*k1, u, p);
   k3 = f(x+h/2*k2, u, p);
   k4 = f(x+h*k3,   u, p);
   x_next = x + h/6*(k1+2*k2+2*k3+k4);
