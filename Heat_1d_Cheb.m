% Heat 1d Cheb - Chebyshev spectral collocation.
%
% Erik Boström, erikbos@kth.se
%   
close all; clear all; clc;
alpha = 0.05;
tau=2.;

N = 30;
[x,DM]=chebdif(N,1);

u0 = exp(-((x-1.1).^2)/0.1)/sqrt(2*pi);

D1 = DM;
D1(1,:) = 0;

dt = 0.00001;
Nt = 10000;

b = zeros(N,1);
b(1) = tau;
u = u0;
for i = 1:Nt
  u =u+dt*DM*(D1*u+b);
  u(end) = alpha;
  plot(x,u,'-'); pause(0.01);
end
