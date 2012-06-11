function prob2 = p_marginalization(M,x0_b,b_table,prob_M,op_fb)
%% Marginalization of prob_M
if nargin < 5
    op_fb = 0;
end

if op_fb ~= 0
    % simultaneous backward and forward
    N = length(M);
    N0_b = length(x0_b);
    prob2 = zeros(2^N0_b,2^N0_b);
    x0_b_M = zeros(N0_b,1);
    two_pow = 2.^(0: N0_b-1);
    two_pow = two_pow';
    for i=1: N0_b
        x0_b_M(i) = find(M==x0_b(i));
    end
    prob_M = reshape(prob_M,[2^N 2^N]);
    for i=1: 2^N
        xp_bs = b_table{i,N};
        xp_bs = xp_bs(x0_b_M);
        xp_i = trans10(xp_bs,two_pow);
        for j=1: 2^N
            xf_bs = b_table{j,N};
            xf_bs = xf_bs(x0_b_M);
            xf_i = trans10(xf_bs,two_pow);
            % fprintf('xp_i=%d xf_i=%d\n',xp_i,xf_i);
            prob2(xp_i,xf_i) = prob2(xp_i,xf_i) + prob_M(i,j);
        end
    end
else
    % backward or forward
    N = length(M);
    N0_b = length(x0_b);
    prob2 = zeros(2^N0_b,1);
    x0_b_M = zeros(N0_b,1);
    two_pow = 2.^(0: N0_b-1);
    two_pow = two_pow';
    for i=1: N0_b
        x0_b_M(i) = find(M==x0_b(i));
    end
    for i=1: 2^N
        x0_bs = b_table{i,N};
        x0_bs = x0_bs(x0_b_M);
        x0_i = trans10(x0_bs,two_pow);
        prob2(x0_i) = prob2(x0_i) + prob_M(i);
    end
    
end