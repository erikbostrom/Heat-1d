% heat 1d FEM
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

alpha = 0.2;
tau=2;

u0 = exp(-((x-1).^2)/0.1)/sqrt(2*pi);

e = -ones(n,1);
K = spdiags([e -2*e e], -1:1, n, n);
K(1,1)=1; K(end,end)=1;
K = K/dx/dx;

e = ones(n,1);
M = spdiags([e 4*e e], -1:1, n, n);
M = (1/6)*M;
M(1,1)=2; M(end,end)=2;
iM = inv(M);

S = zeros(n,1);
S(1) = tau;
S=S/dx;

u = u0';
dt = 0.00001;
Nt = 1000;
I = eye(n);

figure(1);
for t = 0:Nt;
    u = (I - dt*iM*K)*u - dt*iM*S;
    u(end) = alpha;
    plot(x,u); pause(0.01);
end
