function [prob] = sig_prob(x1,C,k,p)

% probability of x1 given x0 (indexed by k)
% x1:the firing pattern of x1, C: subset of x1, k: index of the firing patterns of x0
% p: probability matrix

N = length(C);

prob = 1;
for i=1: N
    j = C(i); % index of x1
    if x1(j) == 1
        prob = prob*p(j,k);
    else
        prob = prob*(1-p(j,k));
    end
end
