function prob = partial_prob_comp(x0_b,x1_b,x1,p,b_table)
% compute the conditional probability
% x1: data, x0_b: partition of x0, x1_b: partition of x1
% p: probability matrix in the whole system

N = log2(size(p,1)); % number of elements in the whole system
two_pow = zeros(N,1);
for i=1: N
    two_pow(i) = 2^(i-1);
end

N0_b = length(x0_b); % number of elements in x0
x0_r = 1: N;
x0_r(x0_b) = []; % the rest of x0
N0_r = length(x0_r); % number of elements in the rest of x0

N1_b = length(x1_b);
x1_r = 1:N;
x1_r(x1_b) = []; % the rest of x1
N1_r = length(x1_r); % number of elements in the rest of x1

prob = zeros(2^N0_b,1); % the coditional entropy p(x0_b(not fixed)|x1_b(fixed))

x1_i1 = sum(two_pow(x1_b).*x1(x1_b));

if N0_r == 0
    x0_i2_vec = 0;
else
    x0_i2_vec = zeros(2^N0_r,1);
    for j=1: 2^N0_r
        % x0_rs = trans2(j-1,N0_r);
        x0_rs = b_table{j,N0_r};
        x0_i2 = sum(two_pow(x0_r).*x0_rs);
        x0_i2_vec(j) = x0_i2;
    end
end

if N1_r == 0
    x1_i2_vec = 0;
else
    x1_i2_vec = zeros(2^N1_r,1);
    for j=1: 2^N1_r
        x1_rs = b_table{j,N1_r};
        x1_i2 = sum(two_pow(x1_r).*x1_rs);
        x1_i2_vec(j) = x1_i2;
    end
end

%% new version
for i=1: 2^N0_b
    % x0_bs = trans2(i-1,N0_b);
    x0_bs = b_table{i,N0_b};
    x0_i1 = sum(two_pow(x0_b).*x0_bs);
    % x0_s(x0_b) = x0_bs;
    for j=1: 2^N0_r % summation over the rest of x0
        % x0_rs = trans2(j-1,N0_r);
        % x0_s(x0_r) = x0_rs;
        % x0_i = trans10(x0_s,two_pow);
        x0_i2 = x0_i2_vec(j);
        x0_i = 1+x0_i1 + x0_i2;
        if N1_r ~= 0
            for k=1: 2^N1_r % summation over the rest of x1
                % x1_rs = trans2(k-1,N1_r);
                % x1_s(x1_r) = x1_rs;
                % x1_i = trans10(x1_s,two_pow);
                x1_i2 = x1_i2_vec(k);
                x1_i = 1+x1_i1+x1_i2;
                % fprintf('x0=%s x1=%s p=%f\n',mat2str(trans2(x0_i-1,N)),mat2str(trans2(x1_i-1,N)),p(x0_i,x1_i));
                prob(i) = prob(i) + p(x0_i,x1_i);
            end
        else
            % x1_i = trans10(x1_s,two_pow);
            x1_i = 1+x1_i1;
            prob(i) = prob(i) + p(x0_i,x1_i);
        end
    end
end

%% old version
% x0_s = zeros(N,1); % the states of x0
% x1_s = zeros(N,1); % the states of x1
% x1_s(x1_b) = x1(x1_b);
% for i=1: 2^N0_b
%     x0_bs = trans2(i-1,N0_b);
%     x0_s(x0_b) = x0_bs;
%     for j=1: 2^N0_r % summation over the rest of x0
%         x0_rs = trans2(j-1,N0_r);
%         x0_s(x0_r) = x0_rs;
%         x0_i = trans10(x0_s,two_pow);
%         if N1_b ~= N
%             for k=1: 2^N1_r % summation over the rest of x1
%                 x1_rs = trans2(k-1,N1_r);
%                 x1_s(x1_r) = x1_rs;
%                 x1_i = trans10(x1_s,two_pow);
%                 prob(i) = prob(i) + p(x0_i,x1_i);
%             end
%         else
%             x1_i = trans10(x1_s,two_pow);
%             prob(i) = prob(i) + p(x0_i,x1_i);
%         end
%     end
% end

%% Normalization
if sum(prob) ~= 0
    prob = prob/sum(prob);
end