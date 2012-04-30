function [T0] = getT0_s0(s0)
	% gets T0 given s0 
	if s0 >= 2/1000
		T0 = -0.03590 - 0.0499*s0 - 0.000112*s0^2;
	else
		T0 = -0.06807*s0;
	end
end