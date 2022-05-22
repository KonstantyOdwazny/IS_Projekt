clear all;clc;
Data = load('cstr.dat');
t = Data(:,1);
u = Data(:,2);
y = Data(:,4);
% plot(t,u,t,y);
% legend('u','y')
% figure()
% plot(t,u)
% title('u')
% figure()
% plot(t,y)
% title('y')
N = 7500;
Tp = 0.1;
n = -N+1:1:N-1;
sigmau = var(u);
avgu = mean(u);
sigmay = var(y);
avgy = mean(y);
[ruu, tauu] = xcorr(u','biased');
[ryy, tauy] = xcorr(y,'biased');

figure()
plot(n*Tp, ryy)