function [muhat, sigmahat, rhohat, Zt, state_index] = calc_markov(T,Tsim,Z,Zprob)
% Anna Cororaton
% simple monte carlo markov chain simulation
% inputs: T - length of data simulated
%         Tsim - number of simulations
%         Z - states
%         Zprob - transition probabilities
% outputs: muhat - mean
%          sigmahat - standard deviation
%          rhohat - first order autocorrelation
% note: can vectorize simulations later

paramhat = zeros(Tsim,3);
display('simulation run...')
for i = 1:Tsim
    
    eps = rand(T,1);
    state_index = zeros(T,1);
    state_index(1,1) = randi(length(Z),1);
    
    for t = 1:T-1
        state_index(t+1,1) = find(cumsum(Zprob(state_index(t,1),:),2)- eps(t+1)>=0,1,'first');
    end
    
    Zt = Z(state_index);
    paramhat(i,1) = nanmean(Zt); % mean 
    paramhat(i,2) = nanstd(Zt);
    Ztl1 = [nan;Zt(1:length(Zt)-1)];
    paramhat(i,3) = nancorr(Zt,Ztl1);
    
    i;
    
end

muhat = nanmean(paramhat(:,1));
sigmahat = nanmean(paramhat(:,2));
rhohat = nanmean(paramhat(:,3));

end