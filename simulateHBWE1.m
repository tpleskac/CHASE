i = 1;
%the gamble
gambles(i).x = 4;
gambles(i).y = 0;
gambles(i).q = .8;
gambles(i).r_x = gambles(i).q;
gambles(i).r_y = 1; 
gambles(i).ev = gambles(i).q.*gambles(i).x + (1-gambles(i).q).*gambles(i).y; 
gambles(i).var = gambles(i).q.*(gambles(i).x-gambles(i).ev )^2 + (1-gambles(i).q).*(gambles(i).y-gambles(i).ev )^2;
gambles(i).sd = gambles(i).var.^.5;

%random certain equivalent
gambles(i).crand = 3;%gambles(i).ev.*.2.*rand + .1.*gambles(i).ev; %a random certaintye equivalent that is equally likely to be below to above

%simulate 1000 players with iCHASE

phiStar = 0; %criterion
pstay = 0; % going to set this at 0 and then move around. 
tau = 1; %standard deviation of start point
thresholdStar = 3; %threshold
alpha = 1; %exponent for utility function (power function)
gamma = 1; %sensitivity for prelec prob. weighting function


pars(1) = phiStar; %criterion
pars(2) = pstay; % going to set this at 0 and then move around. 
pars(3) = tau; %standard deviation of start point
pars(4) = thresholdStar; %threshold
pars(5) = alpha; %exponent for utility function (power function)
pars(6) = gamma; %sensitivity for prelec prob. weighting function

[trialStats, summaryStats] = simiChaseWeightSingle(pars, gambles)