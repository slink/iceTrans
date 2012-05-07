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
etaH = zetaH ./ (3*Pr)^(1/4); % based on equation 28 on page 4
Fp = y(:,2);
phi = y(:,4);
S = y(:,6);
	
figure('Position',[500 300 1.4*400 400]);
plot(etaH,Fp, etaH,phi, etaH,S)
h = xlabel('$\eta, non-dimentional distance downstream$'); set(h, 'interpreter', 'latex');
h1 = legend('F''','$\phi$','S'); set(h1, 'interpreter', 'latex');
set(gcf,'PaperPositionMode','auto');
matlab2tikz('ratioder2Tstore_mesh.tikz',...
 'height', '\figureheight', 'width', '\figurewidth', 'showInfo',false);
	
 %%% need to code the conversion back to u, T, and s for the second plot
 
 