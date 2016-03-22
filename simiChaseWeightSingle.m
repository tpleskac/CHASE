function [trialStats, summaryStats] = simiChaseWeightSingle(pars, gambles)
%pars is a vector of parameters for iCHASE with rank dependent weights on a
%where searching a single gamble and choosing between a sure thing
%gambles is structure with gamble stats
    
%parameters for iCHASE*
phiStar = pars(1); %criterion
pstay = pars(2); % going to set this at 0 and then move around. 
tau = pars(3); %standard deviation of start point
thresholdStar = pars(4); %threshold
alpha = pars(5); %exponent for utility function (power function)
gamma = pars(6); %sensitivity for prelec prob. weighting function
delta = 1; %assuming no shift in weighting function
beta = alpha; %assuming same power function for gains and losses
lambda = 1; %assuming no loss; 
scale = 2;

for sim = 1:10000;
    % simulate a trial from each of the 6 gamble types
    for typeI = 1:length(gambles)
        %set up trial variables
         
        observe = [];%vector storing what observed for trial
        action = []; %vector storing action -1 sample; 0 choose 
        value = []; %vector storing incremental preference for each draw
        pref = []; %vector storing the accumulated preference...
        
        pi_x = prelec(gambles(typeI).r_x,gamma,delta);
        pi_y = (1 - prelec(gambles(typeI).r_x,gamma,delta));
        prospectEU = pi_x.*ptvalue(gambles(typeI).x,alpha,beta,lambda); 
        prospectSD = (pi_x.*(ptvalue(gambles(typeI).x,alpha,beta,lambda) - prospectEU).^2 + ...
                     pi_y.*(ptvalue(gambles(typeI).y,alpha,beta,lambda) - prospectEU).^2).^.5;
        
        %standardizing the criterion and threshold
        phi = phiStar.*scale.*prospectSD; 
        threshold = thresholdStar.*scale.*prospectSD;
        
        %set up sampling weights
        omega_x = (prelec(gambles(typeI).r_x,gamma,delta) - 0)./gambles(typeI).q;
        omega_y = (1 - prelec(gambles(typeI).r_x,gamma,delta))./(1-gambles(typeI).q);
        
        %noise in the starting point. note using a normal distribution
        %because in this model preference is a continuous state, but
        %truncated by threshold
        z = threshold+1; 
        while (z > threshold || z < -threshold)
            z = tau.*randn;
        end;
        pref(1) = z; %assume unbiased preference state. 
        action(1) = -1;
        %the process
        draw = 1; %counter for number of draws.
        while (pref(draw) < threshold && pref(draw) > -threshold); % when pref is between thresholds make a draw, observe, and update preference
            draw = draw + 1;
            %what do I observe
            action(draw) = -1;
            observe(draw) = gambles(typeI).y;%assume going to see outcome y
            omega = omega_y;
            randObserve = rand;
            if (randObserve<gambles(typeI).q)%with probability q see outcome x
                observe(draw) = gambles(typeI).x;
                omega = omega_x;
            end;
            %value observation
            value(draw) = omega.*ptvalue(observe(draw),alpha,beta,lambda)-phi - ptvalue(gambles(typeI).crand,alpha,beta,lambda); 
            %update preference
            randStay = rand;
            if (randStay<pstay)
               pref(draw) = pref(draw-1); 
            else
               pref(draw) = pref(draw-1) + value(draw); 
            end;
        end;

        %the choice
        action(draw) = 0; 
        if (pref(draw) > threshold)
            choice = 1;%upper threshold
        else
            choice = 0;%lower threshold
        end;
        sampleSize = draw - 1; 
        
        trialStats(typeI,sim).observe = observe;%vector storing what observed for trial
        trialStats(typeI,sim).action = action; %vector storing action -1 sample; 0 choose 
        trialStats(typeI,sim).value = value; %vector storing incremental preference for each draw
        trialStats(typeI,sim).pref = pref; %vector storing the accumulated preference...

        trialStats(typeI,sim).choice = choice; 
        trialStats(typeI,sim).sampleSize = sampleSize;
    end;
end;
for typeI = 1:length(gambles)
   %analytical data
    pi_x = prelec(gambles(typeI).r_x,gamma,delta);
    pi_y = (1 - prelec(gambles(typeI).r_x,gamma,delta));
    prospectEU = pi_x.*ptvalue(gambles(typeI).x,alpha,beta,lambda); 
    prospectSD = (pi_x.*(ptvalue(gambles(typeI).x,alpha,beta,lambda) - prospectEU).^2 + ...
                  pi_y.*(ptvalue(gambles(typeI).y,alpha,beta,lambda) - prospectEU).^2).^.5;

    threshold = ceil(thresholdStar); %this needs to be discrete...so round up
    d = (prospectEU - ptvalue(gambles(typeI).crand,alpha,beta,lambda))./(scale.*prospectSD); % note upper threshold is risk (not HEV) 

    [Q, R, Z, I] = CHASEChoiceMatrices(d, pstay, threshold, tau);
    [PrBPrA] = CHASEchoice(Q, R, Z, I);
    [EtBEtA] = CHASEExpDecSample(Q, R, Z, I, PrBPrA);
    summaryStats(typeI).preChoicePr = PrBPrA(2);
    summaryStats(typeI).preMeanSample = PrBPrA*EtBEtA'; 
    
   %simulated data
   theChoice = [trialStats(typeI,:).choice];
   theSampleSize = [trialStats(typeI,:).sampleSize];
   summaryStats(typeI).choiceProp = mean(theChoice);
   sampleStats = [trialStats(typeI,:).sampleSize];
   summaryStats(typeI).medianSample = median(theSampleSize); 
   summaryStats(typeI).meanSample = mean(theSampleSize); 
end;