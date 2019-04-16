close all;
clear;
clc;

%% Damped Spring Example
%
%  dx/dt = v
%  dv/dt = (-k/m)*x + (-c/m)*v + (1/m)*F
%
%  dq/dt = A*q + B*u
%  y = C*q + D*u
%  
%  dq/dt = [dx/dt ; dv/dt]
%  y = x

k = 0.5;
c = 0.1;
m = 1.0;

t = 0:0.1:50;
u = zeros(length(t), 1);

x = [ 1     0    ]';

A = [ 0     1     ;
      -k/m  -c/m ];
   
B = [ 0     ;
      1/m  ];

C = [ 1     0    ];

D = [ 0 ];

sys = ss(A,B,C,D);

%% Kalman Filter initialisation
Q = 1;
R = 1;

rng default
% w = 0.1*randn(length(t),1);
w = zeros(length(t),1);
% v = 0.1*randn(length(t),1);
v = t'./100;

P = B*Q*B';
y = lsim(sys,u+w,t,x);
yv = y+v;

ye = zeros(length(t),1);
ycov = zeros(length(t),1);

plot(t,y);xlabel('time');ylabel('position');

for i = 1:length(t)
  % Measurement update
  Mn = P*C'/(C*P*C'+R);
  x = x + Mn*(yv(i)-C*x);   % x[n|n]
  P = (eye(2)-Mn*C)*P;      % P[n|n]

  ye(i) = C*x;
  errcov(i) = C*P*C';

  % Time update
  x = A*x + B*u(i);        % x[n+1|n]
  P = A*P*A' + B*Q*B';     % P[n+1|n]
end

subplot(211), plot(t,y,'--',t,ye,'-',t,yv,'--')
title('Time-varying Kalman filter response')
xlabel('No. of samples'), ylabel('Output')
subplot(212), plot(t,y-yv,'-.',t,y-ye,'-')
xlabel('No. of samples'), ylabel('Output')