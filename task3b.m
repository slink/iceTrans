clear all
close all
clc

N=4;
v0=NaN(N,1);
v1=NaN(N,1);

for i=1:N
    [k, hil, cp, rhol, rhoi, Pr, Sc, Tinf, T0, s0, sinf, q, rhom0, Tm0, alpha0, g, g1, g2, g3, CGr, u]...
        = getPhysprop(i);
    Le = u{6};
    a =  0.3;
    b =  -2.04;
    
    phiP0bar = b / (Le * (cp * (T0 - Tinf) / hil) * s0 / ((1-s0/1000) * (s0 - sinf)));
    F0 = ((-phiP0bar * cp * (T0 - Tinf))/(hil * (1 - s0 / 1000)));
    
    y0 = [F0, 0, a, 1, phiP0bar, 1, b];
    zeta0 = 0;
    zetaE = 18;
    
    [zetaH, y, anew, bnew] = shootingMethod(zeta0, zetaE, u, a, b, cp, T0, Tinf, hil, s0, sinf);
    
    x=linspace(1e-3,1,100);
    dTdy=mean(y(1,5)*(T0-Tinf)*(3*Pr)^(1/4)./(sqrt(2)*x).*(CGr*x.^3).^(1/4));
    v0(i)=k/((1-s0*1e-3)*rhol*hil)*dTdy;
    v1(i)=v0(i)*rhol*(1-1e-3*s0)/rhoi;
end