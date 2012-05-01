clc, clear all, close all
[k, hil, cp, rhol, Pr, Sc, Tinf, T0, sinf, q, rhom0, Tm0, alpha0, g, g1, g2, g3]... 
	= getPhysprop; 
	
s0 = gets0_T0(T0);

% calculated equations 8 through 11 with s = sinf and T = Tinf
rhom = rhom0*(1+g1*sinf);
alph = alpha0*(1+g2*sinf);
Tm = Tm0*(1+g3*sinf); 

rho = rhom * (1-alph*abs(Tinf-Tm)^q);

R = (Tm - Tinf) / (T0 - Tinf);
A = (g1 * rhom0 * (s0 - sinf)) / rhom;
B = (g2 * alpha0 * (s0 - sinf)) / alph;
Q = (g3 * Tm0 * (s0 - sinf)) / (T0 - Tinf);
P = (g1 * rhom0 * (s0 - sinf)) / (rhom * alph * abs(T0 - Tinf)^q);
Le = Sc / Pr;

% sprintf('-A*100 = %d,\n B*100 = %d,\n -Q =  %d,\n -P = %d,\n R = %d\n', -A*100, B*100, -Q, -P, R)

% initial guess values for a and b
a = 0.123 / (3*Pr)^(-1/4);
b = -5.104;

phiP0bar = b / (Le * (cp * (T0 - Tinf) / hil) * s0 / ((1-s0/1000) * (s0 - sinf)));
F0 = ((-phiP0bar * cp * (T0 - Tinf))/(hil * (1 - s0 / 1000)));

sprintf('phiP0bar = %d,\n F0 = %d,\n', -phiP0bar, F0)

u = {R; A; B; Q; P; Le; q};

y0 = [F0, 0, a, 1, phiP0bar, 1, b];
zeta0 = 0;
zetaE = 5;

% [zetaH, y] = ode23s(@(zetaH, y) ydiff(zetaH, y, u), [zeta0 zetaE], y0);
dxhat = 1e-4; 
jj = 0;
n = 500;
Ns = 20e4;

y = shootingMethodExp(dxhat, y0, Ns, n, u, a, b, cp, T0, Tinf, hil, s0, sinf);

%{
ysave = zeros(Ns/n,7);
y = ydiffExp(dxhat, y0, u);
for ii = 1:Ns
	ynew = ydiffExp(dxhat, y, u);
	if mod(ii,n) == 0
		jj = jj + 1;
		ysave(jj,:) = ynew;
	end
	y = ynew;
end
%}

%{ 
plot(zetaH, y(:,2), zetaH, y(:,6))
[zetaH, y] = shootingMethod(zeta0, zetaE, u, a, b, cp, T0, Tinf, hil, s0, sinf);
%}