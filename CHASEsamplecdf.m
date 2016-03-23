function [PrTBPrTA] = CHASEsamplecdf(Q, R, Z, I, n, PrBPrA, conditional)
%n = round(time./tau);

PrTBPrTA = Z*inv(I-Q)*(I-Q^n)*R./(PrBPrA.^conditional);

i = PrTBPrTA < 0;
PrTBPrTA(i) = 0;
