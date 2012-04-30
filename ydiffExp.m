function ydiff = ydiffExp(dxhat, y, u)
	Ahat = dxhat .* fhat(y, u);
	Bhat = dxhat .* fhat(y + 0.5.*Ahat, u);
	Chat = dxhat .* fhat(y + 0.5.*Bhat, u);
	Dhat = dxhat .* fhat(y + Chat, u);
	
	ydiff = y + (Ahat + 2.*Bhat + 2.*Chat + Dhat)./6;
	
end
