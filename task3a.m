clc, clear all, close all

[k, hil, cp, rhol, rhoi, Pr, Sc, Tinf, T0, s0, sinf, q, rhom0, Tm0, alpha0, g, g1, g2, g3, CGr, u]...
    = getPhysprop(3);
	
Le = u{6};
a =  0.4;
b =  -1.8;
    
phiP0bar = b / (Le * (cp * (T0 - Tinf) / hil) * s0 / ((1-s0/1000) * (s0 - sinf)));
F0 = ((-phiP0bar * cp * (T0 - Tinf))/(hil * (1 - s0 / 1000)));
    
y0 = [F0, 0, a, 1, phiP0bar, 1, b];
zeta0 = 0;
zetaE = 18;
    
[zetaH, y, anew, bnew] = shootingMethod(zeta0, zetaE, u, a, b, cp, T0, Tinf, hil, s0, sinf);
eta = zetaH/(3*Pr)^(1/4);
Fcap = y(:,1);	f = Fcap/(3*Pr)^(3/4);
Fpcap = y(:,2); fp = Fpcap/(3*Pr)^(1/2);
F2pcap = y(:,3); f2p = F2pcap/(3*Pr)^(1/4);
phi = y(:,4);
phip = y(:,5);
S = y(:,6);
Sp = y(:,7);

figure('Position',[500 300 1.4*400 400]);
plot(eta,Fpcap,eta,phi,eta,S)
h = xlabel('$\eta$, non-dimentional distance downstream'); set(h, 'interpreter', 'latex');
ylabel('[N/A]')
h1 = legend('F''','$\phi$','S'); set(h1, 'interpreter', 'latex');
axis([0 max(eta) -0.01 1])
set(gcf,'PaperPositionMode','auto');
matlab2tikz('etaVsF_phi_S_3a.tikz',...
 'height', '\figureheight', 'width', '\figurewidth', 'showInfo',false);
close


%Functions of physical variable
x0 = 1; % y = linspace(0,0.01,100);
yax = sqrt(2)*x0/((CGr*x0^3)^(1/4))*eta;

mu=Pr*k/(cp*1000);
s = sinf+S*(s0-sinf);
T = Tinf+phi*(T0-Tinf);
u = 2/x0*(mu/rhol)*(CGr*x0^3)^(1/2)*fp;

sPlot = s(1e-7<diff(s));
yPlot = yax(1e-7<diff(s));

figure('Position',[500 300 1.4*400 400]);
plot(yPlot * 100,sPlot)
h = ylabel('Salinity, s [${}^o/{}_{oo}$]');  set(h, 'interpreter', 'latex');
xlabel('Distance away from interface, y [cm]')
set(gcf,'PaperPositionMode','auto');
% axis([0 1.1*max(yPlot) 5 11])
% set(gca,'xticklabel',num2str(get(gca,'xtick')'))
% set(gca,'XTickLabel',['0';'0.02';'0.04';'0.06';'0.08';'0.10'])
matlab2tikz('salinity_3a.tikz',...
 'height', '\figureheight', 'width', '\figurewidth', 'showInfo',false);
close

figure('Position',[500 300 1.4*400 400]);
plot(yax,T)
xlim([0 max(yax)])
h = ylabel('T [${}^{\circ}$ C]'); set(h, 'interpreter', 'latex');
xlabel('Distance away from interface, y [m]')
set(gcf,'PaperPositionMode','auto');
matlab2tikz('temperture_3a.tikz',...
 'height', '\figureheight', 'width', '\figurewidth', 'showInfo',false);
close

figure('Position',[500 300 1.4*400 400]);
plot(yax,u*1e3)
xlim([0 max(yax)])
ylabel('u [mm/s]')
xlabel('Distance away from interface, y [m]')
set(gcf,'PaperPositionMode','auto');
matlab2tikz('velocity_3a.tikz',...
 'height', '\figureheight', 'width', '\figurewidth', 'showInfo',false);
close

% Heat transfer coefficient
 
x=linspace(1e-3,1,100);
Grx = CGr*x.^3;
dTdy = y(1,5)*(3*Pr)^(1/4)./(sqrt(2)*x).*(CGr*x.^3).^(1/4);
h = -k*dTdy;

figure('Position',[500 300 1.4*400 400]);
plot(x,h)
xlabel('x, Height of ice [m]')
h = ylabel('h [$W/m^2 K$]'); set(h, 'interpreter', 'latex');
set(gcf,'PaperPositionMode','auto');
matlab2tikz('heat_transfer_coeff_3a.tikz',...
 'height', '\figureheight', 'width', '\figurewidth', 'showInfo',false);
close
 