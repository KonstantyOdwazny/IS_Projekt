clear all;clc;
Data = load('cstr.dat');
t = Data(:,1);
u = Data(:,2);
y = Data(:,4);
N = 7500;
Tp = 0.1;
n = -N+1:1:N-1;
sigmau = var(u);
avgu = mean(u);
sigmay = var(y);
avgy = mean(y);
%% Wykresy
plot(t,u,t,y);
legend('u','y')
figure()
plot(t,u)
title('u')
figure()
plot(t,y)
title('y')
%% Metoda korelacyjna

[ruu, tauu] = xcorr(u','biased');
[ryy, tauy] = xcorr(y,'biased');
[ryu, tauuy] = xcorr(y,u,'biased');
ruu = ruu(7500:end);
ryu = ryu(7500:end);
% figure()
% plot(n*Tp, ryy)

%metoda nieparametryczna
M = 100;

ry = ryu(1:M);
um = u(1:M);
Ruu = zeros(M,M);
for i=1:M
    for j=1:M
        if( i > j)
            Ruu(i,j) = ruu(i-j+1);
        else
            if (j-i ~= 0)
                Ruu(i,j) = ruu(j-i+1);
            else
                Ruu(i,j) = ruu(1);
            end
        end
    end
end
gM = inv(Ruu)*ry*(1/Tp);

hm = [];
tn = 0:1:length(gM)-1;
for i = 1:length(gM)
    hm(i) = Tp*sum(gM(1:i));
end
figure()
stairs(Tp*tn,hm,'r')
%% Metoda analizy widmowej
[ruu, tauu] = xcorr(u,'biased');
[ryu, tauuy] = xcorr(y,u,'biased');
m = 1000;
k = 0:1:m-1;
omegak = 2*pi*k/m;
omega = omegak/Tp;
Gn = fft(y)./fft(u);
G16 = fft(ryu(7500:end))./fft(ruu(7500:end)); 
mod = abs(Gn(1:m));
arg = angle(Gn(1:m))*180/pi;
for i = 1:m
    if arg(i)>2*pi
        arg(i) = arg(i)-360;
    end
end
Lm1 = 20*log10(mod);
figure()
semilogx(omega, Lm1, 'g')
title('Modol')
figure()
semilogx(omega,arg,'k')
title('Angle')

mod16 = abs(G16(1:m));
arg2 = angle(G16(1:m))*180/pi;
for i = 1:m
    if arg2(i)> 10
        arg2(i) = arg2(i)-360;
    end
end
Lm2 = 20*log10(mod16);            
figure()
semilogx(omega, Lm2, 'b')
title('Modol')
figure()
semilogx(omega,arg2,'r')
title('Angle')

%% Metoda LS
                