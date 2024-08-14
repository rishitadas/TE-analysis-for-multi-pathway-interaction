% Computes the entropy of N time series of length T
function [H] = compute_entropy(ts)

T = length(ts); % length of the time series
[l,r,s]=unique(ts,'rows');
p_xs=hist(s,1:length(l))/T; % compute the probability of each combination of realizations
p_xs=p_xs+(p_xs==0); % transform 0 to 1
H=-sum(p_xs.*log2(p_xs)); % compute entropy in bits

end