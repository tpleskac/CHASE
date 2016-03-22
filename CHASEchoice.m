function [PrBPrA] = CHASEchoice(Q, R, Z, I)
PrBPrA = (Z*inv(I-Q)*R);

%PrBPrA = PrBPrA./sum(PrBPrA);
