clc, clear all, close all
[k, hil, cp, rhol, rhoi, Pr, Sc, Tinf, T0, s0, sinf, q, rhom0, Tm0, alpha0, g, g1, g2, g3, CGr, u]...
    = getPhysprop(1);

%sprintf('-A*100 = %d,\n B*100 = %d,\n -Q =  %d,\n -P = %d,\n R = %d\n', -A*100, B*100, -Q, -P, R)

Le = u{6};

% initial guess values for a and b
a =  0.3; % 0.123 / (3*Pr)^(-1/4);
b =  -2.04; % -5.104;

phiP0bar = b / (Le * (cp * (T0 - Tinf) / hil) * s0 / ((1-s0/1000) * (s0 - sinf)));
F0 = ((-phiP0bar * cp * (T0 - Tinf))/(hil * (1 - s0 / 1000)));

%sprintf('phiP0bar = %d,\n F0 = %d,\n', -phiP0bar, F0)

y0 = [F0, 0, a, 1, phiP0bar, 1, b];
zeta0 = 0;
zetaE = 18;

% [zetaH, y] = ode23s(@(zetaH, y) ydiff(zetaH, y, u), [zeta0 zetaE], y0);
% Ns = 2e5;
% dxhat = zetaE / Ns;
% n = Ns / 1000;

%[yExp, anewExp, bnewExp] = shootingMethodExp(dxhat, y0, Ns, n, u, a, b, cp, T0, Tinf, hil, s0, sinf);

[zetaH, y, anew, bnew] = shootingMethod(zeta0, zetaE, u, a, b, cp, T0, Tinf, hil, s0, sinf);

%plot(zetaH,y)
%legend('F','F''','F''''','\phi','\phi''','S','S''')

z=zetaH;
F=y(:,1);
Fp=y(:,2);
F2p=y(:,3);
phi=y(:,4);
phip=y(:,5);
S=y(:,6);
Sp=y(:,7);

figure
plot(z,Fp,z,phi,z,S)
legend('F''','\phi','S')
xlabel('\zeta')

s=sinf+S*(s0-sinf);
T=Tinf+phi*(T0-Tinf);


sPlot=s(1e-7<diff(s));
zPlot=z(1e-7<diff(s));
figure
plot(zPlot,sPlot)
ylabel('s')
%legend('s','T')
xlabel('\zeta')

figure
plot(z,T)
ylabel('T')
%legend('s','T')
xlabel('\zeta')

%Heat transfer coefficient
x=linspace(0,1,100);
mu=Pr*k/cp;
Grx=rhom0/rhom*g*g1*x.^3*(sinf-s0)/(mu/rhol)^2;
h=k*(T0-Tinf)*phip(1)*(3*Pr)^(1/4)./(sqrt(2)*x).*(Grx).^(1/4);

figure
plot(x,h)
xlabel('x')
ylabel('h [W/m^2/K]')



%x=z/(3*Pr)^(1/4)/()




% legend('F','F''','F''''','\phi','\phi''','S','S''')
%{
jj = 0;
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