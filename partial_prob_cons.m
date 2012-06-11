function prob = partial_prob_cons(x0_so,x0_in,x0_out,x1_b,xc,p,b_table,op_fb)

% compute the backward or forward repertoire given x0
% backward: p(xp(not fixed)|x0(fixed)), forward: p(xf(not fixed)|x0(fixed))
% x0_so: fixed state, x0_in: injecting noise in units (maxEnt)
% x0_out: injecting noise in connections (complete noise)
% xc: given current state
% p: transition probability matrix (TPM)
% b_table: table used for converting binary sequences into decimal number
% op_fb: 0 -> forward, 1-> backward

if op_fb == 0
    %% forward
    Nf = length(x1_b);
    if Nf ~= 0
        prob = zeros(2^Nf,1);
        for i=1: 2^Nf
            xf = b_table{i,Nf};
            prob(i) = partial_prob_forward(x0_so,x0_in,x0_out,x1_b,xc,xf,p,b_table);
        end
        prob = prob/sum(prob); % normalization
    else
        prob = [];
    end
else
    %% backward
    Np = length(x0_so);
    if Np ~= 0
        prob = zeros(2^Np,1);
        for i=1: 2^Np
            xp = b_table{i,Np};
            prob(i) = partial_prob_forward(x0_so,x0_in,x0_out,x1_b,xp,xc,p,b_table);
        end
        prob = prob/sum(prob); % normalization
    else
        prob = [];
    end
end