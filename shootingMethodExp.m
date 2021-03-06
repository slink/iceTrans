function [y, anew, bnew] = shootingMethodExp(dxhat, y0, Ns, n, u, a, b, cp, T0, Tinf, hil, s0, sinf)
	% wrap the explicit RK 4 and thus allow for us to then use the shotting method with it
	Le = u{6};
	
	phiP0bar = b / (Le * (cp * (T0 - Tinf) / hil) * s0 / ((1-s0/1000) * (s0 - sinf)));
	F0 = ((-phiP0bar * cp * (T0 - Tinf))/(hil * (1 - s0 / 1000)));
	y0 = [F0, 0, a, 1, phiP0bar, 1, b];
	
	y = ydiffExpRunner(dxhat,y0,u,Ns,n);
	
	% initial errors given the guesses
	Ea = y(end, 3); 
	Eb = y(end, 6);
	
	% initialization for the while loop
	dEaDa = 1;
	dEbDb = 1;
	
	ii = 1;
	% first guess for a new a and b
	anew(ii) = 1.02 * a;
	bnew(ii) = 1.02 * b;
	
	% based on the above, first difference bettween a and anew (also initial values)
	da = 0.02*a;
	db = 0.02*b;
	
	
	while sqrt(dEaDa^2 + dEbDb^2) > 5E-3 % 10^-5

		% ode for anew (1.02 * a)
		phiP0barnewA = b / (Le * (cp * (T0 - Tinf) / hil) * s0 / ((1-s0/1000) * (s0 - sinf)));
		F0newA = ((-phiP0barnewA * cp * (T0 - Tinf))/(hil * (1 - s0 / 1000)));
		y0newA = [F0newA, 0, anew(ii), 1, phiP0barnewA, 1, b];
		ynewA = ydiffExpRunner(dxhat,y0newA,u,Ns,n);
		% [zetaHnewA, ynewA] = ode23s(@(zetaH, y) ydiff(zetaH, y, u), [zeta0 zetaE], y0newA);
		
		% ode for bnew (1.02 * b)
		phiP0barnewB = bnew(ii) / (Le * (cp * (T0 - Tinf) / hil) * s0 / ((1-s0/1000) * (s0 - sinf)));
		F0newB = ((-phiP0barnewB * cp * (T0 - Tinf))/(hil * (1 - s0 / 1000)));
		y0newB = [F0newB, 0, a, 1, phiP0barnewB, 1, bnew(ii)];
		ynewB = ydiffExpRunner(dxhat,y0newB,u,Ns,n);
		% [zetaHnewB, ynewB] = ode23s(@(zetaH, y) ydiff(zetaH, y, u), [zeta0 zetaE], y0newB);
		
		ii = ii + 1
		% if ii == 3
		%	keyboard
		% end
		
		% keyboard
			
		% error for the anew and bnew runs of the ode solver
		EaDa = ynewA(end, 3);
		EbDa = ynewA(end, 6);
		EaDb = ynewB(end, 3);
		EbDb = ynewB(end, 6);
		
		% a and b derivatives based on the anew and bnew runs
		
		dEaDa = (EaDa - Ea) / da;
		dEaDb = (EaDb - Ea) / db;
		dEbDa = (EbDa - Eb) / da;
		dEbDb = (EbDb - Eb) / db;
		
		%{
		dEaDa = (EaDa - Ea) / (0.2*a);
		dEaDb = (EaDb - Ea) / (0.2*b);
		dEbDa = (EbDa - Eb) / (0.2*a);
		dEbDb = (EbDb - Eb) / (0.2*b);
		%}
		
		% new guesses for a and b and their differences
		D = dEaDa*dEbDb - dEaDb*dEbDa;
		da = (-Ea*dEbDb + Eb*dEaDb) / D;
		db = (-Eb*dEaDa + Ea*dEbDa) / D;
		anew(ii) = a + da;
		bnew(ii) = b + db;
		
		% ode run for the old a and b 
		phiP0bar = b / (Le * (cp * (T0 - Tinf) / hil) * s0 / ((1-s0/1000) * (s0 - sinf)));
		F0 = ((-phiP0bar * cp * (T0 - Tinf))/(hil * (1 - s0 / 1000)));
		y0 = [F0, 0, a, 1, phiP0bar, 1, b];
		y = ydiffExpRunner(dxhat,y0,u,Ns,n);
		% [zetaH, y] = ode23s(@(zetaH, y) ydiff(zetaH, y, u), [zeta0 zetaE], y0);
		
		Ea = y(end, 3); 
		Eb = y(end, 6);
		
		a = anew(ii);
		b = bnew(ii);
		
	end	
	
end
