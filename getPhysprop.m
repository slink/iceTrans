function [k, hil, cp, rhol, rhoi, Pr, Sc, Tinf, T0, s0, sinf, q, rhom0, Tm0, alpha0, g, g1, g2, g3, CGr, u]... 
	= getPhysprop(ii)

k = 0.575; % W/m*C
hil = 334.0; % kJ/kg
cp = 4.079; % kJ/kg*C
rhol = 999.97; % kg/m^3
rhoi = 916.7;

Pr = [13.31 13.04 12.79 13.42];
Sc = [2965 2825 2697 3018];
Tinf = [0 1 2 0];
T0 = [-0.478 -0.382 -0.313 -0.875];
sinf = [10 10 10 20];

Tinf = Tinf(ii);
sinf = sinf(ii);
T0 = T0(ii);
Pr = Pr(ii);
Sc = Sc(ii);

q = 1.894816;
rhom0 = 999.9720; % kg/m^3
Tm0 = 4.029325; % C
alpha0 = 9.297173E-6; % C^q
g = 9.807; % m/s^2
g1 = 8.046157E-4;
g2 = -2.839092E-3;
g3 = -5.265509E-2;

%Calculate constants

s0 = gets0_T0(T0);
% calculated equations 8 through 11 with s = sinf and T = Tinf
rhom = rhom0*(1+g1*sinf);
alph = alpha0*(1+g2*sinf);
Tm = Tm0*(1+g3*sinf); 

R = (Tm - Tinf) / (T0 - Tinf);
A = (g1 * rhom0 * (s0 - sinf)) / rhom;
B = (g2 * alpha0 * (s0 - sinf)) / alph;
Q = (g3 * Tm0 * (s0 - sinf)) / (T0 - Tinf);
P = (g1 * rhom0 * (s0 - sinf)) / (rhom * alph * abs(T0 - Tinf)^q);
Le = Sc / Pr;
u = {R; A; B; Q; P; Le; q};

mu=Pr*k/cp;
CGr=rhom0/rhom*g*g1*(sinf-s0)/(mu/rhol)^2;

end
