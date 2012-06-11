function [H H1 H2] = KLD(prob,prob2)
prob2(prob2==0) = 1; % avoid log0 when computing entropy
H1 = - sum(prob.*log2(prob2)) ;

prob(prob==0) = 1;
H2 = - sum(prob.*log2(prob));
H = H1 - H2;