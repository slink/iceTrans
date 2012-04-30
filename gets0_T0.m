function [s0] = gets0_T0(T0)
	% gets s0 given T0 -- inversion of getT0_s0
	if T0 <= -0.1361
		s0 = roots([-0.000112, -0.0499, (-0.03590 - T0)]);
		% root1 = (0.0499 + sqrt(0.0499^2 - 4*0.03590*0.000112))/ (4 * -0.03590);
		% root2 = (0.0499 - sqrt(0.0499^2 - 4*0.03590*0.000112))/ (4 * -0.03590);
		% s0 = [root1 root2];
		s0 = s0(find(s0 > 0));
	else
		s0 = T0 / (-0.06807);
	end
end