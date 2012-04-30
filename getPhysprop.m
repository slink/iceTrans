function [k, hil, cp, rhol, Pr, Sc, Tinf, T0, sinf]...
    =getPhysprop()

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

end
