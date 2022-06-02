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

% [ruu, tauu] = xcorr(u,'biased');
% [ryy, tauy] = xcorr(y,'biased');
% [ryu, tauuy] = xcorr(y,u,'biased');
% ruu = ruu(7500:end);
% ryu = ryu(7500:end);

M = 50;
Ruu = zeros(M,M);
ryu = zeros(M,1);
tn = 0:1:M-1;
for i=1:M
    for j=1:M
        Ruu(i,j) = Covar([u,u],j-i);
    end
    ryu(i,1) = Covar([y,u],i-1); 
end
gM = inv(Ruu)*ryu;
figure()
plot(tn*Tp,gM);

%% Metoda analizy widmowej
[ruu, tauu] = xcorr(u,'unbiased');
[ryu, tauuy] = xcorr(y,u,'unbiased');
Mw = N/10;
k = 0:1:N-1;
Np = (N/2)-1;
omegak = 2*pi*k/N; % [rad]
omega = omegak/Tp; %[rad/s]
wh = zeros(N,1);
% Okno Hamminga
for i=0:N-1
    if (i < Mw)
        wh(i+1,1) = 0.5*(1 + cos(i*pi/Mw));
    end
end
ru = zeros(N,1);
ry = zeros(N,1);
for i=0:N-1
    ru(i+1,1) = Covar([u,u],i)*wh(i+1);
    ry(i+1,1) = Covar([y,u],i)*wh(i+1);
end

Gn = fft(y)./fft(u);
Gw = fft(ryu(7500:end))./fft(ruu(7500:end));
mod = abs(Gn(1:Np));
arg = angle(Gn(1:Np))*180/pi;
for i = 1:Np
    if arg(i)>2*pi
        arg(i) = arg(i)-360;
    end
end
mod1 = abs(Gw(1:Np));
arg1 = angle(Gw(1:Np))*180/pi;
for i = 1:Np
    if arg1(i)>2*pi
        arg1(i) = arg1(i)-360;
    end
end
Lm = 20*log10(mod);
Lm1 = 20*log10(mod1);
figure()
semilogx(omega(1:Np), Lm, 'b')
title('Modol')
figure()
semilogx(omega(1:Np),arg,'k')
title('Angle')




%% Metoda LS
ilosc_parametrow = 2;
yest = y(1:4000);
uest = u(1:4000);
ywer = y(4000:end);
uwer = u(4000:end);
fi = zeros(length(uest), ilosc_parametrow);
for i = ilosc_parametrow:length(yest)
%     fi(i,1) = -1*yest(i-1);
%     fi(i,2) = -1*yest(i-2);
%     fi(i,3) = uest(i-1);
    for j = 1:ilosc_parametrow-1
        fi(i,j) = -1*yest(i-j);
    end
    fi(i,ilosc_parametrow) = uest(i-1);
end
p = inv(fi' * fi)*fi' * yest;
ym = zeros(1,length(ywer));
yp = zeros(1,length(ywer));
for i = ilosc_parametrow:length(ywer)
%     ym(i) = -1*p(1)*ym(i-1) - p(2)*ym(i-2) + p(3)*uwer(i-1);
%     yp(i) = -1*p(1)*ywer(i-1)- p(2)*ywer(i-2) + p(3)*uwer(i-1);
    for j=1:ilosc_parametrow-1
        ym(i) = ym(i)- p(j)*ym(i-j);
        yp(i) = yp(i)- p(j)*ywer(i-j);
    end
    ym(i) = ym(i)+ p(ilosc_parametrow)*uwer(i-1);
    yp(i) = yp(i)+ p(ilosc_parametrow)*uwer(i-1);
end
hold on
plot(y)
plot(yp,'g')
title('LS')
legend('ywer - zmierzone','yp - odpowiedz predykowana')
hold off
%% Zmienne instrumentalne
yest = y(1:4000);
uest = u(1:4000);
x =zeros(1,length(yest));
z = zeros(length(yest),2);
sumx=0;
for i=1:length(yest)
    if i==1
        x(i) = p(1)*yest(i)+p(2)*uest(i);
        z(i,1) = (z((i),1) - x(i));
        z(i,2) = (z((i),2) - uest(i));
    else
        x(i) = p(1)*yest(i-1)+p(2)*uest(i-1);
        z(i,1) = (z((i-1),1) - x(i-1));
        z(i,2) = (z((i-1),2) - uest(i-1));
    end
end
piv = (inv(z' * fi))*z' * yest;
ymiv = zeros(1,length(ywer));
ypiv = zeros(1,length(ywer));
for i = 2:length(ywer)
    ymiv(i) = -1*piv(1)*ymiv(i-1) + piv(2)*uwer(i-1);
    ypiv(i) = -1*piv(1)*ywer(i-1) + piv(2)*uwer(i-1);
end
figure()
hold on
plot(ywer)
plot(ymiv,'k')
title('IV')
legend('ywer - zmierzone','yiv - odpowiedz predykowana')
hold off