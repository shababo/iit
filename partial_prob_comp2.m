function prob = partial_prob_comp2(x0_so,x0_in,x0_out,x1_b,source,p,b_table,op_fb)

if op_fb == 0 || op_fb == 2
    % forward computation
    x0 = source(x0_so);
    N1_b = length(x1_b);
    prob = zeros(2^N1_b,1);
    for i=1: 2^N1_b % states of target x1
        x1 = b_table{i,N1_b};
        prob(i) = partial_prob_forward(x0_so,x0_in,x0_out,x1_b,x0,x1,p,b_table);
    end
elseif op_fb == 1 || op_fb == 3
    % backward computation
    x1 = source(x1_b);
    N0_b = length(x0_so);
    prob = zeros(2^N0_b,1);
    for i=1: 2^N0_b % state of target x0
        x0 = b_table{i,N0_b};
        prob(i) = partial_prob_forward(x0_so,x0_in,x0_out,x1_b,x0,x1,p,b_table);
    end
end

% Normalization (in forward, probably not needed. just in case)
if op_fb ~= 3
    if sum(prob) ~= 0
        prob = prob/sum(prob);
    end
end