function [dydxhat] = fhat(y, u)
	
	% a 4th order RK 
	R = u{1}; A = u{2};
	B = u{3}; Q = u{4};
	P = u{5}; Le = u{6};
	q = u{7};
	% C3 = 1 !!!
	ydiff(1) = y(2);
	ydiff(2) = y(3);
	ydiff(3) = 1 / P * ((1 + A * y(6)) * (1 + B * y(6))...
	 			* abs(y(4) - R - Q*y(6))^q - abs(R)^q - P*y(6));
	ydiff(4) = y(5);
	ydiff(5) = -y(1) * y(5);
	ydiff(6) = y(7);
	ydiff(7) = -Le * y(1) * y(7);
	
	dydxhat = ydiff;
	
	
end