function [EtBEtA] = CHASEExpDecSample(Q, R, Z, I, PrBPrA)

EtBEtA = (Z*(I-Q)^(-2)*R)./PrBPrA;