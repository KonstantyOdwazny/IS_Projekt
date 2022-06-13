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
% plot(t,u,t,y);
% legend('u','y')
% figure()
% plot(t,u)
% title('u')
% figure()
% plot(t,y)
% title('y')
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
title('Metoda korelacyjna');

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
    if arg(i) > 2*pi
        arg(i) = arg(i)-360;
    end
end
mod1 = abs(Gw(1:Np));
arg1 = angle(Gw(1:Np))*180/pi;
for i = 1:Np
    arg(i) = arg(i)+180;
      
 
end

Lm = 20*log10(mod);
Lm1 = 20*log10(mod1);
figure()
subplot(2,1,1)
semilogx(omega(1:Np), Lm, 'b')
title('Modol')
ylabel('20*log(mod)')
subplot(2,1,2)
semilogx(omega(1:Np),arg,'k')
title('Angle')
ylabel('Deg')




%% Metoda LS
ilosc_parametrow = 2;
yest = y(1:4000);
uest = u(1:4000);
ywer = y(4000:end);
uwer = u(4000:end);
% Model statyczny
fis = zeros(length(uest), 2);
for i = 1:length(yest)
    fis(i,1) = uest(i);
    fis(i,2) = 1;
end
ps = pinv(fis) * yest;
yms = zeros(1,length(ywer));
for i = 1:length(ywer)
    yms(i) = ps(1)*uwer(i) + ps(2);
end
figure()
hold on
plot(ywer)
plot(yms, 'r')
title('Model statyczny')
legend('y - zmierzone','ym - odpowiedz symulatora')
hold off
% Model FIR
fir = zeros(length(uest), 1);
for i = 2:length(yest)
    fir(i,1) = uest(i-1);
%     fir(i,2) = uest(i-2);
end
pr = pinv(fir) * yest;
ymr = zeros(1,length(ywer));
ymr(1) = 440;
for i = 2:length(ywer)
    ymr(i) = pr(1)*uwer(i-1); %+ pr(2)*uwer(i-2);
end
% figure()
% hold on
% plot(ywer)
% plot(ymr, 'r')
% title('Model FIR')
% legend('y - zmierzone','ym - odpowiedz symulatora')
% hold off
% Model oscylacyjny
fi = zeros(length(uest), 3);
for i = 3:length(yest)
    fi(i,1) = -1*yest(i-1);
    fi(i,2) = -1*yest(i-2);
    fi(i,3) = uest(i-2);
end
p = pinv(fi) * yest;
ym = zeros(1,length(ywer));
yp = zeros(1,length(ywer));
ym(1) = 440;
yp(1) = 440;
ym(2) = 440;
yp(2) = 440;
for i = 3:length(ywer)
    ym(i) = -1*p(1)*ym(i-1) - p(2)*ym(i-2) + p(3)*uwer(i-2);
    yp(i) = -1*p(1)*ywer(i-1)- p(2)*ywer(i-2) + p(3)*uwer(i-2);
end
% figure()
% hold on
% plot(ywer)
% plot(yp,'g')
% plot(ym, 'r')
% title('Model oscylacyjny')
% legend('y - zmierzone','yp - odpowiedz predykowana','ym - odpowiedz symulatora')
% hold off
% ARX 2 parametry plus offset
fi2 = zeros(length(uest), 3);
for i = 2:length(yest)
    fi2(i,1) = -1*yest(i-1);
    fi2(i,2) = uest(i-1);
    fi2(i,3) = 1;
end
p2 = pinv(fi2) * yest;
ym2 = zeros(1,length(ywer));
yp2 = zeros(1,length(ywer));
ym2(1) = 440;
yp2(1) = 440;
for i = 2:length(ywer)
    ym2(i) = -1*p2(1)*ym2(i-1) + p2(2)*uwer(i-1) + p2(3);
    yp2(i) = -1*p2(1)*ywer(i-1) + p2(2)*uwer(i-1) + p2(3);
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
fi3 = zeros(length(uest), 5);
for i = 3:length(yest)
    fi3(i,1) = -1*yest(i-1);
    fi3(i,2) = -1*yest(i-2);
    fi3(i,3) = uest(i-1);
    fi3(i,4) = uest(i-2);
    fi3(i,5) = 1;
end
p3 = pinv(fi3) * yest;
ym3 = zeros(1,length(ywer));
yp3 = zeros(1,length(ywer));
ym3(1) = 440;
yp3(1) = 440;
ym3(2) = 440;
yp3(2) = 440;
for i = 3:length(ywer)
    ym3(i) = -1*p3(1)*ym3(i-1) - p3(2)*ym3(i-2) + p3(3)*uwer(i-1) + p3(4)*uwer(i-2) + p3(5);
    yp3(i) = -1*p3(1)*ywer(i-1) - p3(2)*ywer(i-2) + p3(3)*uwer(i-1) + p3(4)*uwer(i-2) + p3(5);
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
fi4 = zeros(length(uest), 7);
for i = 4:length(yest)
    fi4(i,1) = -1*yest(i-1);
    fi4(i,2) = -1*yest(i-2);
    fi4(i,3) = -1*yest(i-3);
    fi4(i,4) = uest(i-1);
    fi4(i,5) = uest(i-2);
    fi4(i,6) = uest(i-3);
    fi4(i,7) = 1;
end
p4 =pinv(fi4) * yest;
ym4 = zeros(1,length(ywer));
yp4 = zeros(1,length(ywer));
ym4(1) = 440;
yp4(1) = 440;
ym4(2) = 440;
yp4(2) = 440;
ym4(3) = 440;
yp4(3) = 440;
for i = 4:length(ywer)
    ym4(i) = -1*p4(1)*ym4(i-1) - p4(2)*ym4(i-2)-p4(3)*ym4(i-3) + p4(4)*uwer(i-1) + p4(5)*uwer(i-2) + p4(6)*uwer(i-3) + p4(7);
    yp4(i) = -1*p4(1)*ywer(i-1) - p4(2)*ywer(i-2)- p4(3)*ywer(i-3) + p4(4)*uwer(i-1) + p4(5)*uwer(i-2)+ p4(6)*uwer(i-3) + p4(7);
end
figure()
hold on
plot(ywer)
plot(yp4,'g')
plot(ym4, 'r')
title('Model ARX 6 parametrow metoda LS')
legend('y - zmierzone','yp - odpowiedz predykowana','ym - odpowiedz symulatora')
hold off


%% Bledy

eps = ywer' - yms;
Vs = (1/length(ywer)) * eps*eps';

epr = ywer' - ymr;
Vr = (1/length(ywer)) * epr*epr';

Vp = zeros(1,3);
Vm = zeros(1,3);
% ep1 = ywer' - yp;
% Vp(1) = (1/length(ywer)) * ep1*ep1';
% em1 = ywer' - ym;
% Vm(1) = (1/length(ywer)) * em1*em1';
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
tm = [2 4 6];
figure()
plot(tm,log(Vp),tm,log(Vm))
legend('Vp','Vm');
xlabel('m')
ylabel('ln(V)')
%% Zmienne instrumentalne
yest = y(1:4000);
uest = u(1:4000);
ywer = y(4000:end);
uwer = u(4000:end);
x = zeros(1,length(yest));
z = zeros(length(yest),3);
sumx=0;
% IV 2 parametry
for i=2:length(yest)
    x(i) = -1*p2(1)*x(i-1)+p2(2)*uest(i-1);
    z(i,1) = -1*x(i-1);
    z(i,2) = uest(i-1); 
    z(i,3) = 1;
end
piv = (inv(z' * fi2))*z' * yest;
ymiv = zeros(1,length(ywer));
ypiv = zeros(1,length(ywer));
ymiv(1) = 440;
ypiv(1) = 440;
for i = 2:length(ywer)
    ymiv(i) = -1*piv(1)*ymiv(i-1) + piv(2)*uwer(i-1) + piv(3);
    ypiv(i) = -1*piv(1)*ywer(i-1) + piv(2)*uwer(i-1) + piv(3);
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
z2 = zeros(length(yest),5);
for i=3:length(yest)
    x2(i) = -1*p3(1)*x2(i-1) - p3(2)*x2(i-2)+ p3(3)*uest(i-1) + p3(4)*uest(i-2);
    z2(i,1) = -1*x2(i-1);
    z2(i,2) = -1*x2(i-2);
    z2(i,3) = uest(i-1); 
    z2(i,4) = uest(i-2);
    z2(i,5) = 1;
end
piv2 = (inv(z2' * fi3))*z2' * yest;
ymiv2 = zeros(1,length(ywer));
ypiv2 = zeros(1,length(ywer));
ymiv2(1) = 440;
ypiv2(1) = 440;
ymiv2(2) = 440;
ypiv2(2) = 440;
for i = 3:length(ywer)
    ymiv2(i) = -1*piv2(1)*ymiv2(i-1) - piv2(2)*ymiv2(i-2) + piv2(3)*uwer(i-1) + piv2(4)*uwer(i-2) + piv2(5);
    ypiv2(i) = -1*piv2(1)*ywer(i-1) - piv2(2)*ywer(i-2)  + piv2(3)*uwer(i-1) + piv2(4)*uwer(i-2) + piv2(5);
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
z3 = zeros(length(yest),7);
for i=4:length(yest)
    x3(i) = -1*p4(1)*x3(i-1) - p4(2)*x3(i-2) - p4(3)*x3(i-3) + p4(4)*uest(i-1) + p4(5)*uest(i-2) + p4(6)*uest(i-3);
    z3(i,1) = -1*x3(i-1);
    z3(i,2) = -1*x3(i-2);
    z3(i,3) = -1*x3(i-3);
    z3(i,4) = uest(i-1); 
    z3(i,5) = uest(i-2); 
    z3(i,6) = uest(i-3); 
    z3(i,7) = 1;
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
    ymiv3(i) = -1*piv3(1)*ymiv3(i-1) - piv3(2)*ymiv3(i-2) - piv3(3)*ymiv3(i-3) + piv3(4)*uwer(i-1) + piv3(5)*uwer(i-2) + piv3(6)*uwer(i-3) + piv3(7);
    ypiv3(i) = -1*piv3(1)*ywer(i-1) - piv3(2)*ywer(i-2)- piv3(3)*ywer(i-3)  + piv3(4)*uwer(i-1) + piv3(5)*uwer(i-2)+ piv3(6)*uwer(i-3) + piv3(7);
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
% Gld = tf([-0.0088],[1,-1],Tp);
% Gld2 = tf([-0.121,0.1221],[1,-1.9252,0.9254],Tp);
% Gld3 = tf([-0.119,0.1715,-0.0521],[1,-2.67,2.43362,-0.7628],Tp);
% 
% Gid = tf([0.0231],[1,-0.9947],Tp);
% Gid2 = tf([-0.1214,0.1217],[1,-1.8989,0.8990],Tp);
% Gid3 = tf([-0.1196,0.1661,-0.0465],[1,-2.6252,2.3410,-0.7158],Tp);
% 
% r1 = roots([1,-1]);
% r2 = roots([1,-1.9252,0.9254]);
% r3 = roots([1,-2.67,2.43362,-0.7628]);
% zer1 = 0;
% zer2 = roots([-0.121,0.1221]);
% zer3 = roots([-0.119,0.1715,-0.0521]);
% 
% riv1 = roots([1,-0.9947]);
% riv2 = roots([1,-1.8989,0.8990]);
% riv3 = roots([1,-2.6252,2.3410,-0.7158]);
% zeriv1 = 0;
% zeriv2 = roots([-0.1214,0.1217]);
% zeriv3 = roots([-0.1196,0.1661,-0.0465]);
%% Porownanie elastycznosci
% figure()
% hold on
% plot(tn*Tp,gM,'g');
% impulse(Gld)
% impulse(Gld2)
% impulse(Gld3)
% impulse(Gid)
% impulse(Gid2)
% impulse(Gid3)
% axis([0,3,-4,5])
% hold off
% title('Porownanie metody korelacyjnej z metoda LS')
% legend('metoda analizy korealcyjnej','Gls - 2 parametry','Gls2 - 4 parametry','Gls3 - 6 parametrow','Giv - 2 parametry','Giv2 - 4 parametry','Giv3 - 6 parametrow')
% 
% 
% figure()
% hold on
% semilogx(omega(1:Np), Lm, 'k')
% bode(Gld)
% bode(Gld2)
% bode(Gld3)
% bode(Gid)
% bode(Gid2)
% bode(Gid3)
% semilogx(omega(1:Np),arg,'k')
% hold off
% title('Porownanie charakterystyk Bodego')
% legend('Gls - 2 parametry','Gls2 - 4 parametry','Gls3 - 6 parametrow','Giv - 2 parametry','Giv2 - 4 parametry','Giv3 - 6 parametrow','Metoda analizy widmowej')

figure()
plot(tmiv,log(Vmiv),tmiv,log(Vm))
legend('Vmiv','Vmls');
xlabel('m')
ylabel('ln(V)')
re = zeros(1,length(em));
re2 = zeros(1,length(em));
re3 = zeros(1,length(em));
lin = (1.96/sqrt(length(em)))*ones(1,length(em));
lin2 = (-1.96/sqrt(length(em)))*ones(1,length(em));
for i = 1:length(em)
    re(i) = Covar([em',em'],i-1)/Covar([em',em'],0);
    re2(i) = Covar([em2',em2'],i-1)/Covar([em2',em2'],0);
    re3(i) = Covar([em3',em3'],i-1)/Covar([em3',em3'],0);
end
figure()
hold on
plot(re,'b')
plot(re2,'r')
plot(re3,'g')
plot(lin,'--k')
plot(lin2,'--k')
hold off
title('Test bialosci bledow')
legend('r - LS 2 parametry','r1 - LS 4 parametry','r2 - LS 6 parametrow','linia bledu')


s = zeros(1,length(em));
s2 = zeros(1,length(em));
s3 = zeros(1,length(em));
for i = 1:length(em)
    s(i) = Covar([em',uwer],i-1)/sqrt(Covar([em',em'],0)*Covar([uwer,uwer],0));
    s2(i) = Covar([em2',uwer],i-1)/sqrt(Covar([em2',em2'],0)*Covar([uwer,uwer],0));
    s3(i) = Covar([em3',uwer],i-1)/sqrt(Covar([em3',em3'],0)*Covar([uwer,uwer],0));
end
figure()
hold on
plot(s,'b')
plot(s2,'r')
plot(s3,'g')
plot(lin,'--k')
plot(lin2,'--k')
hold off
title('Test skolerowania bledow resztowych')
legend('s - LS 2 parametry','s1 - LS 4 parametry','s2 - LS 6 parametrow','linia bledu')

%% Porownanie zlozonosci modelu
dp = length(p2);
dp2 = length(p3);
dp3 = length(p4);
dpiv = length(piv);
dpiv2 = length(piv2);
dpiv3 = length(piv3);
FPE = zeros(2,3);
AIC = zeros(2,3);
SIC = zeros(2,3);

FPE(1,1) = Vm(1)*(1+dp/N)/(1 - dp/N);
FPE(1,2) = Vm(2)*(1+dp2/N)/(1 - dp2/N);
FPE(1,3) = Vm(3)*(1+dp3/N)/(1 - dp3/N);
FPE(2,1) = Vmiv(1)*(1+dpiv/N)/(1 - dpiv/N);
FPE(2,2) = Vmiv(2)*(1+dpiv2/N)/(1 - dpiv2/N);
FPE(2,3) = Vmiv(3)*(1+dpiv3/N)/(1 - dpiv3/N);

AIC(1,1) = N*log(Vm(1)) + 2*dp;
AIC(1,2) = N*log(Vm(2)) + 2*dp2;
AIC(1,3) = N*log(Vm(3)) + 2*dp3;
AIC(2,1) = N*log(Vmiv(1)) + 2*dpiv;
AIC(2,2) = N*log(Vmiv(2)) + 2*dpiv2;
AIC(2,3) = N*log(Vmiv(3)) + 2*dpiv3;

SIC(1,1) = N*log(Vm(1)) + 2*dp*log(N);
SIC(1,2) = N*log(Vm(2)) + 2*dp2*log(N);
SIC(1,3) = N*log(Vm(3)) + 2*dp3*log(N);
SIC(2,1) = N*log(Vmiv(1)) + 2*dpiv*log(N);
SIC(2,2) = N*log(Vmiv(2)) + 2*dpiv2*log(N);
SIC(2,3) = N*log(Vmiv(3)) + 2*dpiv3*log(N);

Mi = (fi2' * fi2)/4000;
Mi2 = (fi3' * fi3)/4000;
Mi3 = (fi4' * fi4)/4000;
Miv = (z' * fi2)/4000;
Miv2 = (z2' * fi3)/4000;
Miv3 = (z3' * fi4)/4000;

condm1 = sqrt(max(eig(Mi)))/sqrt(min(eig(Mi)));
condm2 = sqrt(max(eig(Mi2)))/sqrt(min(eig(Mi2)));
condm3 = sqrt(max(eig(Mi3)))/sqrt(min(eig(Mi3)));
condmiv = sqrt(max(eig(Miv)))/sqrt(min(eig(Miv)));
condmiv2 = sqrt(max(eig(Miv2)))/sqrt(min(eig(Miv2)));
condmiv3 = sqrt(max(eig(Miv3)))/sqrt(min(eig(Miv3)));

es = zeros(length(yest),1);
sigma = 0;
es1 = zeros(length(yest),1);
sigma1 = 0;
es2 = zeros(length(yest),1);
sigma2 = 0;
es3 = zeros(length(yest),1);
sigma3 = 0;
es4 = zeros(length(yest),1);
sigma4 = 0;
es5 = zeros(length(yest),1);
sigma5 = 0;
for i=1:length(yest)
    es(i) = yest(i) - fi2(i,:) * p2;
    es1(i) = yest(i) - fi3(i,:) * p3;
    es2(i) = yest(i) - fi4(i,:) * p4;
    es3(i) = yest(i) - fi2(i,:) * piv;
    es4(i) = yest(i) - fi3(i,:) * piv2;
    es5(i) = yest(i) - fi4(i,:) * piv3;
end
sigma = sum(es.^2)/(length(yest)-4);
cov = sigma * inv(fi2' * fi2);
sigma1 = sum(es1.^2)/(length(yest)-4);
cov1 = sigma1 * inv(fi3' * fi3);
sigma2 = sum(es2.^2)/(length(yest)-4);
cov2 = sigma2 * inv(fi4' * fi4);
sigma3 = sum(es3.^2)/(length(yest)-4);
cov3 = sigma3 * inv(z' * fi2);
sigma4 = sum(es4.^2)/(length(yest)-4);
cov4 = sigma4 * inv(z2' * fi3);
sigma5 = sum(es5.^2)/(length(yest)-4);
cov5 = sigma5 * inv(z3' * fi4);

% Dodac po jednej kolumnie

pu = [p2(1) - 1.96*sqrt(cov(1,1)),p2(1) + 1.96*sqrt(cov(1,1));
    p2(2) - 1.96*sqrt(cov(2,2)),p2(2) + 1.96*sqrt(cov(2,2))
     p2(3) - 1.96*sqrt(cov(3,3)),p2(3) + 1.96*sqrt(cov(3,3))];
pu2 = [p3(1) - 1.96*sqrt(cov1(1,1)),p3(1) + 1.96*sqrt(cov1(1,1));
    p3(2) - 1.96*sqrt(cov1(2,2)),p3(2) + 1.96*sqrt(cov1(2,2));
    p3(3) - 1.96*sqrt(cov1(3,3)),p3(3) + 1.96*sqrt(cov1(3,3));
    p3(4) - 1.96*sqrt(cov1(4,4)),p3(4) + 1.96*sqrt(cov1(4,4))
    p3(5) - 1.96*sqrt(cov1(5,5)),p3(5) + 1.96*sqrt(cov1(5,5))];
pu3 = [p4(1) - 1.96*sqrt(cov2(1,1)),p4(1) + 1.96*sqrt(cov2(1,1));
    p4(2) - 1.96*sqrt(cov2(2,2)),p4(2) + 1.96*sqrt(cov2(2,2));
    p4(3) - 1.96*sqrt(cov2(3,3)),p4(3) + 1.96*sqrt(cov2(3,3));
    p4(4) - 1.96*sqrt(cov2(4,4)),p4(4) + 1.96*sqrt(cov2(4,4));
    p4(5) - 1.96*sqrt(cov2(5,5)),p4(5) + 1.96*sqrt(cov2(5,5));
    p4(6) - 1.96*sqrt(cov2(6,6)),p4(6) + 1.96*sqrt(cov2(6,6))
     p4(7) - 1.96*sqrt(cov2(7,7)),p4(7) + 1.96*sqrt(cov2(7,7))];

puiv = [piv(1) - 1.96*sqrt(cov3(1,1)),piv(1) + 1.96*sqrt(cov3(1,1));
    piv(2) - 1.96*sqrt(cov3(2,2)),piv(2) + 1.96*sqrt(cov3(2,2))
    piv(3) - 1.96*sqrt(cov3(3,3)),piv(3) + 1.96*sqrt(cov3(3,3))];
puiv2 = [piv2(1) - 1.96*sqrt(cov4(1,1)),piv2(1) + 1.96*sqrt(cov4(1,1));
    piv2(2) - 1.96*sqrt(cov4(2,2)),piv2(2) + 1.96*sqrt(cov4(2,2));
    piv2(3) - 1.96*sqrt(cov4(3,3)),piv2(3) + 1.96*sqrt(cov4(3,3));
    piv2(4) - 1.96*sqrt(cov4(4,4)),piv2(4) + 1.96*sqrt(cov4(4,4))
    piv2(5) - 1.96*sqrt(cov4(5,5)),piv2(5) + 1.96*sqrt(cov4(5,5))];
puiv3 = [piv3(1) - 1.96*sqrt(cov5(1,1)),piv3(1) + 1.96*sqrt(cov5(1,1));
    piv3(2) - 1.96*sqrt(cov5(2,2)),piv3(2) + 1.96*sqrt(cov5(2,2));
    piv3(3) - 1.96*sqrt(cov5(3,3)),piv3(3) + 1.96*sqrt(cov5(3,3));
    piv3(4) - 1.96*sqrt(cov5(4,4)),piv3(4) + 1.96*sqrt(cov5(4,4));
    piv3(5) - 1.96*sqrt(cov5(5,5)),piv3(5) + 1.96*sqrt(cov5(5,5));
    piv3(6) - 1.96*sqrt(cov5(6,6)),piv3(6) + 1.96*sqrt(cov5(6,6))
    piv3(7) - 1.96*sqrt(cov5(7,7)),piv3(7) + 1.96*sqrt(cov5(7,7))];


