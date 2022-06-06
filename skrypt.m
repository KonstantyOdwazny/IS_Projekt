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
% Model inercji pierwszego rzedu
fi = zeros(length(uest), 2);
for i = 2:length(yest)
    fi(i,1) = -1*yest(i-1);
    fi(i,2) = 1;
end
p = pinv(fi) * yest;
ym = zeros(1,length(ywer));
yp = zeros(1,length(ywer));
ym(1) = 440;
yp(1) = 440;
for i = 2:length(ywer)
    ym(i) = -1*p(1)*ym(i-1) + p(2)*uwer(i);
    yp(i) = -1*p(1)*ywer(i-1) + p(2)*uwer(i);
end
figure()
hold on
plot(ywer)
plot(yp,'g')
plot(ym, 'r')
title('Model inercji 1 rzedu metoda LS')
legend('y - zmierzone','yp - odpowiedz predykowana','ym - odpowiedz symulatora')
hold off
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
Vp = zeros(1,4);
Vm = zeros(1,4);
ep1 = ywer' - yp;
Vp(1) = (1/length(ywer)) * ep1*ep1';
em1 = ywer' - ym;
Vm(1) = (1/length(ywer)) * em1*em1';
ep = ywer' - yp2;
Vp(2) = (1/length(ywer)) * ep*ep';
em = ywer' - ym2;
Vm(2) = (1/length(ywer)) * em*em';
ep2 = ywer' - yp3;
Vp(3) = (1/length(ywer)) * ep2*ep2';
em2 = ywer' - ym3;
Vm(3) = (1/length(ywer)) * em2*em2';
ep3 = ywer' - yp4;
Vp(4) = (1/length(ywer)) * ep3*ep3';
em3 = ywer' - ym4;
Vm(4) = (1/length(ywer)) * em3*em3';
tm = [1 2 4 6];
figure()
plot(tm,Vp,tm,Vm)
legend('Vp','Vm');
xlabel('m')
ylabel('V')
%% Zmienne instrumentalne
yest = y(1:4000);
uest = u(1:4000);
ywer = y(4000:end);
uwer = u(4000:end);
x = zeros(1,length(yest));
z = zeros(length(yest),2);
sumx=0;
% IV 2 parametry
for i=2:length(yest)
    x(i) = -1*p2(1)*x(i-1)+p2(2)*uest(i-1);
    z(i,1) = -1*x(i-1);
    z(i,2) = uest(i-1); 
end
piv = (inv(z' * fi2))*z' * yest;
ymiv = zeros(1,length(ywer));
ypiv = zeros(1,length(ywer));
ymiv(1) = 440;
ypiv(1) = 440;
for i = 2:length(ywer)
    ymiv(i) = -1*piv(1)*ymiv(i-1) + piv(2)*uwer(i-1);
    ypiv(i) = -1*piv(1)*ywer(i-1) + piv(2)*uwer(i-1);
end
figure()
hold on
plot(ywer)
plot(ymiv,'k')
plot(ypiv)
title('IV 2 parametry')
legend('ywer - zmierzone','ymiv - odpowiedz symulatora','ypiv - odpowiedz predyktora')
hold off
% IV 4 parametry
x2 = zeros(1,length(yest));
z2 = zeros(length(yest),4);
for i=3:length(yest)
    x2(i) = -1*p3(1)*x2(i-1) - p3(2)*x2(i-2)+ p3(3)*uest(i-1) + p3(4)*uest(i-2);
    z2(i,1) = -1*x2(i-1);
    z2(i,2) = -1*x2(i-2);
    z2(i,3) = uest(i-1); 
    z2(i,4) = uest(i-2); 
end
piv2 = (inv(z2' * fi3))*z2' * yest;
ymiv2 = zeros(1,length(ywer));
ypiv2 = zeros(1,length(ywer));
ymiv2(1) = 440;
ypiv2(1) = 440;
ymiv2(2) = 440;
ypiv2(2) = 440;
for i = 3:length(ywer)
    ymiv2(i) = -1*piv2(1)*ymiv2(i-1) - piv2(2)*ymiv2(i-2) + piv2(3)*uwer(i-1) + piv2(4)*uwer(i-2);
    ypiv2(i) = -1*piv2(1)*ywer(i-1) - piv2(2)*ywer(i-2)  + piv2(3)*uwer(i-1) + piv2(4)*uwer(i-2);
end
figure()
hold on
plot(ywer)
plot(ymiv2,'k')
plot(ypiv2)
title('IV 4 parametry')
legend('ywer - zmierzone','ymiv - odpowiedz symulatora','ypiv - odpowiedz predyktora')
hold off
% IV 6 parametry
x3 = zeros(1,length(yest));
z3 = zeros(length(yest),6);
for i=4:length(yest)
    x3(i) = -1*p4(1)*x3(i-1) - p4(2)*x3(i-2) - p4(3)*x3(i-3) + p4(4)*uest(i-1) + p4(5)*uest(i-2) + p4(6)*uest(i-3);
    z3(i,1) = -1*x3(i-1);
    z3(i,2) = -1*x3(i-2);
    z3(i,3) = -1*x3(i-3);
    z3(i,4) = uest(i-1); 
    z3(i,5) = uest(i-2); 
    z3(i,6) = uest(i-3); 
end
piv3 = (inv(z3' * fi4))*z3' * yest;
ymiv3 = zeros(1,length(ywer));
ypiv3 = zeros(1,length(ywer));
ymiv3(1) = 440;
ypiv3(1) = 440;
ymiv3(2) = 440;
ypiv3(2) = 440;
ymiv3(3) = 440;
ypiv3(3) = 440;
for i = 4:length(ywer)
    ymiv3(i) = -1*piv3(1)*ymiv3(i-1) - piv3(2)*ymiv3(i-2) - piv3(3)*ymiv3(i-3) + piv3(4)*uwer(i-1) + piv3(5)*uwer(i-2) + piv3(6)*uwer(i-3);
    ypiv3(i) = -1*piv3(1)*ywer(i-1) - piv3(2)*ywer(i-2)- piv3(3)*ywer(i-3)  + piv3(4)*uwer(i-1) + piv3(5)*uwer(i-2)+ piv3(6)*uwer(i-3);
end
figure()
hold on
plot(ywer)
plot(ymiv3,'k')
plot(ypiv3)
title('IV 6 parametry')
legend('ywer - zmierzone','ymiv - odpowiedz symulatora','ypiv - odpowiedz predyktora')
hold off

Vpiv = zeros(1,3);
Vmiv = zeros(1,3);
epiv1 = ywer' - ypiv;
Vpiv(1) = (1/length(ywer)) * epiv1*epiv1';
emiv1 = ywer' - ymiv;
Vmiv(1) = (1/length(ywer)) * emiv1*emiv1';
epiv = ywer' - ypiv2;
Vpiv(2) = (1/length(ywer)) * epiv*epiv';
emiv = ywer' - ymiv2;
Vmiv(2) = (1/length(ywer)) * emiv*emiv';
epiv2 = ywer' - ypiv3;
Vpiv(3) = (1/length(ywer)) * epiv2*epiv2';
emiv2 = ywer' - ymiv3;
Vmiv(3) = (1/length(ywer)) * emiv2*emiv2';
tmiv = [2 4 6];
figure()
plot(tmiv,log(Vpiv),tmiv,log(Vmiv))
legend('Vpiv','Vmiv');
xlabel('m')
ylabel('ln(V)')
%% Transmitancje
Gld = tf([-0.0088],[1,-1],Tp);
Gld2 = tf([-0.121,0.1221],[1,-1.9252,0.9254],Tp);
Gld3 = tf([-0.119,0.1715,-0.0521],[1,-2.67,2.43362,-0.7628],Tp);

Gid = tf([0.0231],[1,-0.9947],Tp);
Gid2 = tf([-0.1214,0.1217],[1,-1.8989,0.8990],Tp);
Gid3 = tf([-0.1196,0.1661,-0.0465],[1,-2.6252,2.3410,-0.7158],Tp);