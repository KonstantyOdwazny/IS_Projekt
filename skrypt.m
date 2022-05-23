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