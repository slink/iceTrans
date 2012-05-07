clc, clear all, close all

[k, hil, cp, rhol, rhoi, Pr, Sc, Tinf, T0, s0, sinf, q, rhom0, Tm0, alpha0, g, g1, g2, g3, CGr, u]...
    = getPhysprop(3);
	
Le = u{6};
a =  0.3;
b =  -2.04;
    
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
y = sqrt(2)*x0/(CGr*x0^3)^(1/4)*eta;

s = sinf+S*(s0-sinf);
T = Tinf+phi*(T0-Tinf);
mu = Pr*k/cp;
u = 2/x0*(mu/rhol)*(CGr*x0^3)^(1/2)*fp;

sPlot = s(1e-7<diff(s));
yPlot = y(1e-7<diff(s));

figure('Position',[500 300 1.4*400 400]);
plot(yPlot,sPlot)
h = ylabel('Salinity, s [${}^o/{}_{oo}$]');  set(h, 'interpreter', 'latex');
xlabel('Distance away from interface, y [m]')
set(gcf,'PaperPositionMode','auto');
set(gca,'xticklabel',num2str(get(gca,'xtick')'))
% fix based on comment at   ---- tried ' scale only axis'
% http://win.ua.ac.be/~nschloe/comment/255#comment-255
% 'x tick label style', '{/pgf/number format/fixed}',...
% 'xticklabel style', '{/pgf/number format/precision=3}',...
% '/pgfplots/try min ticks', 1,
%'extraAxisOptions', ,...
matlab2tikz('salinity_3a.tikz',...
 'height', '\figureheight', 'width', '\figurewidth', 'showInfo',false);


figure('Position',[500 300 1.4*400 400]);
plot(y,T)
h = ylabel('T [${}^{\circ}$ C]'); set(h, 'interpreter', 'latex');
xlabel('Distance away from interface, y [m]')
set(gcf,'PaperPositionMode','auto');
matlab2tikz('temperture_3a.tikz',...
 'height', '\figureheight', 'width', '\figurewidth', 'showInfo',false);
close

figure('Position',[500 300 1.4*400 400]);
plot(y,u*1e3)
ylabel('u [mm/s]')
xlabel('Distance away from interface, y [m]')
set(gcf,'PaperPositionMode','auto');
matlab2tikz('velocity_3a.tikz',...
 'height', '\figureheight', 'width', '\figurewidth', 'showInfo',false);
close

% Heat transfer coefficient
x = linspace(0,1,100);
mu = Pr*k/cp;
Grx = CGr*x.^3;
h = k*(T0-Tinf)*phip(1)*(3*Pr)^(1/4)./(sqrt(2)*x).*(Grx).^(1/4);

figure('Position',[500 300 1.4*400 400]);
plot(x,h)
xlabel('x, Height of ice [m]')
h = ylabel('h [$W/m^2 K$]'); set(h, 'interpreter', 'latex');
set(gcf,'PaperPositionMode','auto');
matlab2tikz('heat_transfer_coeff_3a.tikz',...
 'height', '\figureheight', 'width', '\figurewidth', 'showInfo',false);
close
 
 