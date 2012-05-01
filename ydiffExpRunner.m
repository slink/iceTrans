function ysave = ydiffExpRunner(dxhat,y0,u,Ns,n)
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
end