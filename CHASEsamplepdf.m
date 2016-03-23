function [PrTBPrTA] = CHASEsamplepdf(Q, R, Z,n, PrBPrA)
%n = round(time./tau);
PrTBPrTA = Z*Q^(n-1)*R./(PrBPrA);
