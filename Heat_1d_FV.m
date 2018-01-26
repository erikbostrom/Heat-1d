% heat 1d FV
% u_t = u_y(tau), tau=u_y, 0<x<1, t>0,
% tau(0,t)=taumodel, u(1,t)=alpha
%
% Erik Boström, erikbos@kth.se
%   
close all; clear all; clc;

n = 30;
lx = 1;
dx = lx/(n-1);
x = 0:dx:lx;

% Boundary coditions
taumodel = 0;
alpha = 0.05;

% Initial condition
u0 = exp(-(x.^2)/0.1)/sqrt(2*pi);
u0 = u0';

e = ones(n,1);
A = spdiags([e -2*e e], -1:1, n, n);
A(1,1) = -1; A(1,2)= 1;
A(end,end) = -3; A(end,end-1)= 1;
A = A/dx/dx;

S = zeros(n,1);
S(1) = taumodel/dx;
S(end) = 2*alpha/dx/dx;

u = u0;
dt = 0.0001;
Nt = 1000;

% Time stepping with Euler explicit.
figure(1);
for t = 0:Nt;
    u = u + dt*A*u + dt*S;
    plot(x,u); pause(0.01);
end
