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
   x_next = x + h*k2;
