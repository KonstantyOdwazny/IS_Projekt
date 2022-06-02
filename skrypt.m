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
% % Model inercji pierwszego rzedu
% fi = zeros(length(uest), 2);
% for i = 2:length(yest)
%     fi(i,1) = -1*yest(i-1);
%     fi(i,2) = uest(i);
% end
% p = pinv(fi) * yest;
% ym = zeros(1,length(ywer));
% yp = zeros(1,length(ywer));
% ym(1) = 440;
% yp(1) = 440;
% for i = 2:length(ywer)
%     ym(i) = -1*p(1)*ym(i-1) + p(2)*uwer(i);
%     yp(i) = -1*p(1)*ywer(i-1) + p(2)*uwer(i);
% end
% figure()
% hold on
% plot(ywer)
% plot(yp,'g')
% plot(ym, 'r')
% title('Model inercji 1 rzedu metoda LS')
% legend('y - zmierzone','yp - odpowiedz predykowana','ym - odpowiedz symulatora')
% hold off
% ARX 2 parametry
fi2 = zeros(length(uest), 2);
for i = 2:length(yest)
    fi2(i,1) = -1*yest(i-1);
    fi2(i,2) = uest(i-1);
end
p2 = pinv(fi2) * yest;
ym2 = zeros(1,length(ywer));
yp2 = zeros(1,length(ywer));
ym2(1) = 440;
yp2(1) = 440;
for i = 2:length(ywer)
    ym2(i) = -1*p2(1)*ym2(i-1) + p2(2)*uwer(i-1);
    yp2(i) = -1*p2(1)*ywer(i-1) + p2(2)*uwer(i-1);
end
figure()
hold on
plot(ywer)
plot(yp2,'g')
plot(ym2, 'r')
title('Model ARX 2 parametry metoda LS')
legend('y - zmierzone','yp - odpowiedz predykowana','ym - odpowiedz symulatora')
hold off
% ARX 4 parametry
fi3 = zeros(length(uest), 4);
for i = 3:length(yest)
    fi3(i,1) = -1*yest(i-1);
    fi3(i,2) = -1*yest(i-2);
    fi3(i,3) = uest(i-1);
    fi3(i,4) = uest(i-2);
end
p3 = inv(fi3' * fi3)*fi3' * yest;
ym3 = zeros(1,length(ywer));
yp3 = zeros(1,length(ywer));
ym3(1) = 440;
yp3(1) = 440;
ym3(2) = 440;
yp3(2) = 440;
for i = 3:length(ywer)
    ym3(i) = -1*p3(1)*ym3(i-1) - p3(2)*ym3(i-2) + p3(3)*uwer(i-1) + p3(4)*uwer(i-2);
    yp3(i) = -1*p3(1)*ywer(i-1) - p3(2)*ywer(i-2) + p3(3)*uwer(i-1) + p3(4)*uwer(i-2);
end
figure()
hold on
plot(ywer)
plot(yp3,'g')
plot(ym3, 'r')
title('Model ARX 4 parametry metoda LS')
legend('y - zmierzone','yp - odpowiedz predykowana','ym - odpowiedz symulatora')
hold off
% ARX 6 parametry
fi4 = zeros(length(uest), 6);
for i = 4:length(yest)
    fi4(i,1) = -1*yest(i-1);
    fi4(i,2) = -1*yest(i-2);
    fi4(i,3) = -1*yest(i-3);
    fi4(i,4) = uest(i-1);
    fi4(i,5) = uest(i-2);
    fi4(i,6) = uest(i-3);
end
p4 = inv(fi4' * fi4)*fi4' * yest;
ym4 = zeros(1,length(ywer));
yp4 = zeros(1,length(ywer));
ym4(1) = 440;
yp4(1) = 440;
ym4(2) = 440;
yp4(2) = 440;
ym4(3) = 440;
yp4(3) = 440;
for i = 4:length(ywer)
    ym4(i) = -1*p4(1)*ym4(i-1) - p4(2)*ym4(i-2)-p4(3)*ym4(i-3) + p4(4)*uwer(i-1) + p4(5)*uwer(i-2) + p4(6)*uwer(i-3);
    yp4(i) = -1*p4(1)*ywer(i-1) - p4(2)*ywer(i-2)- p4(3)*ywer(i-3) + p4(4)*uwer(i-1) + p4(5)*uwer(i-2)+ p4(6)*uwer(i-3);
end
figure()
hold on
plot(ywer)
plot(yp4,'g')
plot(ym4, 'r')
title('Model ARX 6 parametrow metoda LS')
legend('y - zmierzone','yp - odpowiedz predykowana','ym - odpowiedz symulatora')
hold off
Vp = zeros(1,3);
Vm = zeros(1,3);
ep = ywer' - yp2;
Vp(1) = (1/length(ywer)) * ep*ep';
em = ywer' - ym2;
Vm(1) = (1/length(ywer)) * em*em';
ep2 = ywer' - yp3;
Vp(2) = (1/length(ywer)) * ep2*ep2';
em2 = ywer' - ym3;
Vm(2) = (1/length(ywer)) * em2*em2';
ep3 = ywer' - yp4;
Vp(3) = (1/length(ywer)) * ep3*ep3';
em3 = ywer' - ym4;
Vm(3) = (1/length(ywer)) * em3*em3';
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