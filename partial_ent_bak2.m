function [H prob] = partial_ent_bak2(x0_b,x1_b,x1,p,C_j)

% compute the conditional entropy given a partition of x0 and x1
% x1: data, C_x1: data_index, x0_b: partition of x0, x1_b: partition of x1,
% p: probability matrix

if nargin > 4
    op = 1; % K-L divergence
else
    op = 0; % Difference of entropy
    C_j = [];
end

N = size(p,1);

N0_b = length(x0_b); % number of elements in x0
x0_r = 1: N;
x0_r(x0_b) = []; % the rest of x0
N0_r = length(x0_r);

prob_temp = zeros(2^N,1);
prob2_temp = zeros(2^N,1);

for k=1: 2^N
    j = rem(k,2^N0_r);
    if j==0
        j = 2^N0_r;
    end
    i = (k-j)/2^N0_r + 1;
    
    x0_bst = trans2(i-1,N0_b);
    % x0_bst
    x0 = zeros(N,1);
    x0(x0_b) = x0_bst;
    
    x0_rst = trans2(j-1,N0_r);
    % x0_rst
    x0(x0_r) = x0_rst; % rest
    l = trans10(x0);
    % temp = 1;
    %         for k=1: N_d
    %             if x1(x1_b(k)) == 1
    %                 temp = temp*p(x1_b(k),l);
    %             else
    %                 temp = temp*(1-p(x1_b(k),l));
    %             end
    %         end
    prob_temp(k) = sig_prob(x1,x1_b,l,p);
    % pause;
    % x1_b(k)
    
    if op == 1
        prob2_temp(k) = sig_prob(x1,C_j,l,p);
    end
    % fprintf('i=%d j=%d temp=%f ',i,j,temp);
    % pause;
end

prob = zeros(2^N0_b,1);
prob2 = zeros(2^N0_b,1);
for i=1: 2^N0_b
    j = (i-1)*2^N0_r+1:i*2^N0_r;
    prob(i) = sum(prob_temp(j));
    prob2(i) = sum(prob2_temp(j));
end

prob = prob/sum(prob);
prob(prob==0) = 1;

if op == 0
    H = -sum(prob.*log2(prob)); % conditional entropy H(x0_1|x1_1)
else
    prob2 = prob2/sum(prob2);
    prob2(prob2==0) = 1;
    H = -sum(prob2.*log2(prob));
end

display(prob)
display(prob2)
pause;