function [H prob] = partial_ent(x0_b,x1_b,x1,p,C_j)

% x1: data, C_x1: data_index, x0_b: partition of x0, x1_b: partition of x1,
% p: probability matrix

if nargin > 4
    op = 1; % K-L divergence
else
    op = 0; % Difference of entropy
end

N = size(p,1);

N0_b = length(x0_b); % number of elements in x0
x0_r = 1: N;
x0_r(x0_b) = [];
N0_r = length(x0_r);

prob = zeros(2^N0_b,1);
prob2 = zeros(2^N0_b,1);

for i=1: 2^N0_b
    x0_bst = trans2(i-1,N0_b);
    % x0_bst
    x0 = zeros(N,1);
    x0(x0_b) = x0_bst;
    
    for j=1: 2^N0_r
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
        temp = sig_prob(x1,x1_b,l,p);
        % pause;
        % x1_b(k)
        prob(i) = prob(i) + temp;
        
        if op == 1
            temp2 = sig_prob(x1,C_j,l,p);
            prob2(i) = prob2(i) + temp2;
        end
        % fprintf('i=%d j=%d temp=%f ',i,j,temp);
        % pause;
    end    
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
