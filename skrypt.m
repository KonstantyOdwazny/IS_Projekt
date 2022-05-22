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
ruu = ruu(7500:end);
% figure()
% plot(n*Tp, ryy)

%metoda nieparametryczna
M = 10;
k = 1:1:2000;
um = u(1:M);
Ruu = zeros(M,M);
for i=1:M
    for j=1:M
        if( i > j)
            Ruu(i,j) = ruu(i);
        else
            if (j-i ~= 0)
                Ruu(i,j) = ruu(j-i);
            else
                Ruu(i,j) = ruu(1);
            end
        end
    end
end

            
               
            
                