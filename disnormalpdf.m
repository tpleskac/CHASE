function [p] = disnormalpdf(mu, stdev, x)
xlow = x-.5; 
xup = x+.5; 
qlow = normcdf(xlow,mu,stdev);
qup = normcdf(xup,mu,stdev);
p = qup - qlow;
p = p./(normcdf(x(length(x))+.5,mu,stdev)-normcdf(x(1)-.5,mu,stdev));
