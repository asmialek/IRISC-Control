close all;
clear;
clc;

A = [1.1269   -0.4940    0.1129;
     1.0000         0         0;
          0    1.0000         0];

B = [-0.3832;
      0.5919;
      0.5191];

C = [1 0 0];

t = [0:100]';
u = sin(t/5);

Q = 1; 
R = 1;

t = [0:100]';
u = sin(t/5);

n = length(t);
rng default
w = sqrt(Q)*randn(n,1);
v = sqrt(R)*randn(n,1);
% w = zeros(n,1);
% v = t./50;

sys = ss(A,B,C,0,-1);
y = lsim(sys,u+w);      
yv = y + v;

P = B*Q*B';         % Initial error covariance
x = zeros(3,1);     % Initial condition on the state
ye = zeros(length(t),1);
ycov = zeros(length(t),1); 

for i = 1:length(t)
  % Measurement update
  Mn = P*C'/(C*P*C'+R)
  x = x + Mn*(yv(i)-C*x);   % x[n|n]
  P = (eye(3)-Mn*C)*P;      % P[n|n]

  ye(i) = C*x;
  errcov(i) = C*P*C';

  % Time update
  x = A*x + B*u(i);        % x[n+1|n]
  P = A*P*A' + B*Q*B';     % P[n+1|n]
end

Plant = ss(A,[B B],C,0,-1,'inputname',{'u' 'w'},'outputname','y');
[kalmf,L,P,M] = kalman(Plant,Q,R);

subplot(211), plot(t,y,'--',t,ye,'-',t,yv,'--')
title('Time-varying Kalman filter response')
xlabel('No. of samples'), ylabel('Output')
subplot(212), plot(t,y-yv,'-.',t,y-ye,'-')
xlabel('No. of samples'), ylabel('Output')