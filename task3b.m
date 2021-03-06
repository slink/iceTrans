clear all
close all
clc

N=4;
v0=NaN(N,1);
v1=NaN(N,1);
hvec=NaN(N,1);
Tinfvec=NaN(N,1);
CGrvec=NaN(N,1);

%Temperature

for i=1:N
    [k, hil, cp, rhol, rhoi, Pr, Sc, Tinf, T0, s0, sinf, q, rhom0, Tm0, alpha0, g, g1, g2, g3, CGr, u]...
        = getPhysprop(i);
    Le = u{6};
    a =  0.3;
    b =  -2.04;
    Tinfvec(i)=Tinf;
    
    phiP0bar = b / (Le * (cp * (T0 - Tinf) / hil) * s0 / ((1-s0/1000) * (s0 - sinf)));
    F0 = ((-phiP0bar * cp * (T0 - Tinf))/(hil * (1 - s0 / 1000)));
    
    y0 = [F0, 0, a, 1, phiP0bar, 1, b];
    zeta0 = 0;
    zetaE = 18;
    
    [zetaH, y, anew, bnew] = shootingMethod(zeta0, zetaE, u, a, b, cp, T0, Tinf, hil, s0, sinf);
    
    x=linspace(1e-3,1,100);
    dTdy=mean(y(1,5)*(3*Pr)^(1/4)./(sqrt(2)*x).*(CGr*x.^3).^(1/4)); % took out (T0-Tinf)
    hvec(i)= -k*dTdy;
    %v0(i)=k/((1-s0*1e-3)*rhol*hil*1000)*abs(dTdy); % since dTdy is negtive and hil is in KJ
    %v1(i)=v0(i)*rhol/rhoi; % took out (1-1e-3*s0)
    v1(i)=-k*abs(dTdy)/rhoi/(hil*1000)*(T0-Tinf);
    CGrvec(i)=CGr;
end

v1=v1*3600*24*100;

figure('Position',[500 300 1.4*400 400]);
plot(Tinfvec(1:3),hvec(1:3),'-s')
h1 = xlabel('$T_\infty [^{\circ}C]$ '); set(h1, 'interpreter', 'latex');
h2 = ylabel('$h [W/m^2/s]$'); set(h2, 'interpreter', 'latex');
%h1 = legend('F''','$\phi$','S'); set(h1, 'interpreter', 'latex');
xlim([0 2])
set(gcf,'PaperPositionMode','auto');
matlab2tikz('h_Tinf_3b.tikz',...
 'height', '\figureheight', 'width', '\figurewidth', 'showInfo',false);
%close

figure('Position',[500 300 1.4*400 400]);
plot(Tinfvec(1:3),v1(1:3),'-s')
h1 = xlabel('$T_\infty$ '); set(h1, 'interpreter', 'latex');
h2 = ylabel('$\dot{M}$ [cm/day]'); set(h2, 'interpreter', 'latex');
%h1 = legend('F''','$\phi$','S'); set(h1, 'interpreter', 'latex');
xlim([0 2])
set(gcf,'PaperPositionMode','auto');
matlab2tikz('Mdot_Tinf_3b.tikz',...
 'height', '\figureheight', 'width', '\figurewidth', 'showInfo',false);
%close


figure('Position',[500 300 1.4*400 400]);
plot([10 20],hvec([1 4]),'-s')
h1 = xlabel('$s_\infty$ '); set(h1, 'interpreter', 'latex');
h2 = ylabel('$h [W/m^2/s]$'); set(h2, 'interpreter', 'latex');
set(gcf,'PaperPositionMode','auto');
matlab2tikz('h_sinf_3b.tikz',...
 'height', '\figureheight', 'width', '\figurewidth', 'showInfo',false);
%close

figure('Position',[500 300 1.4*400 400]);
plot([10 20],v1([1 4]),'-s')
h1 = xlabel('$s_\infty$ '); set(h1, 'interpreter', 'latex');
h2 = ylabel('$\dot{M}$ [cm/day]'); set(h2, 'interpreter', 'latex');
set(gcf,'PaperPositionMode','auto');
matlab2tikz('Mdot_sinf_3b.tikz',...
 'height', '\figureheight', 'width', '\figurewidth', 'showInfo',false);
%close

%Change in melt rate from T_inf=0 to T_inf=2
dv1T=(v1(3)-v1(1))/2;

%Change in melt rate from s_inf=10 to s_inf=20
dv1s=(v1(4)-v1(1))/10;


%Part 3c, compare Grashoff number with pg. 372 in Keys. Turbulence for Gr >
%10^9-10^10
CGrvec
