function [zetaH, y, anew, bnew] = shootingMethod(zeta0, zetaE, u, a, b, cp, T0, Tinf, hil, s0, sinf)

Le = u{6};
opt=odeset('RelTol',1e-5,'AbsTol',1e-8);

% initialization for the while loop
%dEaDa = 1;
%dEbDb = 1;

% based on the above, first difference bettween a and anew (also initial values)
da = 0.02*a;
db = 0.02*b;

% first guess for a new a and b
anew = a+da;
bnew = b+db;

% original values for the a and b passed to the shooting method
phiP0bar = b / (Le * (cp * (T0 - Tinf) / hil) * s0 / ((1-s0/1000) * (s0 - sinf)));
F0 = ((-phiP0bar * cp * (T0 - Tinf))/(hil * (1 - s0 / 1000)));
y0 = [F0, 0, a, 1, phiP0bar, 1, b];
[zetaH, y] = ode23s(@(zetaH, y) ydiff(zetaH, y, u), [zeta0 zetaE], y0,opt);

% initial errors given the guesses
Ea = y(end, 3);
Eb = y(end, 6);

% ode for anew - uses perturbed value of a and old value for b
phiP0barnewA = b / (Le * (cp * (T0 - Tinf) / hil) * s0 / ((1-s0/1000) * (s0 - sinf)));
F0newA = ((-phiP0barnewA * cp * (T0 - Tinf))/(hil * (1 - s0 / 1000)));
y0newA = [F0newA, 0, anew, 1, phiP0barnewA, 1, b];
[~, ynewA] = ode23s(@(zetaH, y) ydiff(zetaH, y, u), [zeta0 zetaE], y0newA,opt);

% ode for bnew - uses perturbed value of b and old value for a
phiP0barnewB = bnew / (Le * (cp * (T0 - Tinf) / hil) * s0 / ((1-s0/1000) * (s0 - sinf)));
F0newB = ((-phiP0barnewB * cp * (T0 - Tinf))/(hil * (1 - s0 / 1000)));
y0newB = [F0newB, 0, a, 1, phiP0barnewB, 1, bnew];
[~, ynewB] = ode23s(@(zetaH, y) ydiff(zetaH, y, u), [zeta0 zetaE], y0newB,opt);

a=anew;
b=bnew;

while sqrt(Ea^2 + Eb^2) > 1e-5
    sqrt(Ea^2 + Eb^2)
    % error for the anew and bnew runs of the ode solver
    EaDa = ynewA(end, 3);
    EaDb = ynewB(end, 3);
    
    EbDa = ynewA(end, 6);
    EbDb = ynewB(end, 6);
    
    % a and b derivatives based on the anew and bnew runs
    dEaDa = (EaDa - Ea) / da;
    dEaDb = (EaDb - Ea) / db;
    dEbDa = (EbDa - Eb) / da;
    dEbDb = (EbDb - Eb) / db;
    
    % new guesses for a and b and their differences
    D = dEaDa*dEbDb - dEaDb*dEbDa;
    da = (-Ea*dEbDb + Eb*dEaDb) / D;
    db = (-Eb*dEaDa + Ea*dEbDa) / D;
    
    %Over relax the parameters
    da=da*0.3;
    db=db*0.3;
    anew = a + da;
    bnew = b + db;
    
    % ode run for the old a and b
    phiP0bar = b / (Le * (cp * (T0 - Tinf) / hil) * s0 / ((1-s0/1000) * (s0 - sinf)));
    F0 = ((-phiP0bar * cp * (T0 - Tinf))/(hil * (1 - s0 / 1000)));
    y0 = [F0, 0, a, 1, phiP0bar, 1, b];
    [zetaH, y] = ode23s(@(zetaH, y) ydiff(zetaH, y, u), [zeta0 zetaE], y0,opt);
    
    Ea = y(end, 3);
    Eb = y(end, 6);
    
    % ode for anew - uses perturbed value of a and old value for b
    phiP0barnewA = b / (Le * (cp * (T0 - Tinf) / hil) * s0 / ((1-s0/1000) * (s0 - sinf)));
    F0newA = ((-phiP0barnewA * cp * (T0 - Tinf))/(hil * (1 - s0 / 1000)));
    y0newA = [F0newA, 0, anew, 1, phiP0barnewA, 1, b];
    [~, ynewA] = ode23s(@(zetaH, y) ydiff(zetaH, y, u), [zeta0 zetaE], y0newA,opt);
    
    % ode for bnew - uses perturbed value of b and old value for a
    phiP0barnewB = bnew / (Le * (cp * (T0 - Tinf) / hil) * s0 / ((1-s0/1000) * (s0 - sinf)));
    F0newB = ((-phiP0barnewB * cp * (T0 - Tinf))/(hil * (1 - s0 / 1000)));
    y0newB = [F0newB, 0, a, 1, phiP0barnewB, 1, bnew];
    [~, ynewB] = ode23s(@(zetaH, y) ydiff(zetaH, y, u), [zeta0 zetaE], y0newB,opt);
    
    a=anew;
    b=bnew;
end

end