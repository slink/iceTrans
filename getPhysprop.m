function [k, hil, cp, rhol, Pr, Sc, Tinf, T0, sinf, q, rhom0, Tm0, alpha0, g, g1, g2, g3]... 
	= getPhysprop()

k = 0.575; % W/m*C
hil = 334.0; % kJ/kg
cp = 4.079; % kJ/kg*C
rhol = 999.97; % kg/m^3

Pr = [13.31 13.04 12.79 13.42];
Sc = [2965 2825 2697 3018];
Tinf = [0 1 2 0];
T0 = [-0.478 -0.382 -0.313 -0.875];
sinf = [10 10 10 20];

ii = 1;
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

end
