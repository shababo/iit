function [H_max p] = forward_H(M,p,b_table,op_fb)

%% computing the entropy over potential repertoire of x1

N = length(M);

if op_fb == 0
    p_x1 = zeros(2^N,1);
    for j=1: 2^N
        x0 = trans2(j-1,N);
        % maxent
        % p_x1 = p_x1 + partial_prob_comp(M,M,x0,p,b_table,2,M);
        % complete noise
        p_x1 = p_x1 + partial_prob_comp([],M,x0,p,b_table,0,M,M);
    end
    p_x1 = p_x1/2^N;
    
    p = p_x1;
    
    p_x1(p_x1==0) = 1;
    H_max = -sum(p_x1.*log2(p_x1));
    
elseif op_fb == 3
    N_max = size(p,2);
    T = 1;
    p_x1 = zeros(2^N,2^N);
    
    % maxent
%     x0_so = M;
%     x0_in = [];
%     x0_out = pick_rest(1:N_max,M);
    
    % complete noise
    x0_so = [];
    x0_in = [];
    x0_out = 1:N;
    
    x1_in = M;
    x1_out = pick_rest(1:N_max,M);
    for j=1: 2^N_max
        source = trans2(j-1,N_max);
        p_c = patial_prob_comp_bf(T,x0_so,x0_in,x0_out,x1_in,x1_out,source,p,b_table);
        p_x1 = p_x1 + p_c;
    end
    Norm = 2^N_max;
    
    p_x1 = p_x1/Norm;
    p_x1 = p_x1(:);
    
    p = p_x1;
    
    % fprintf('%s\n',mat2str(p));
    % pause;
    
    p_x1(p_x1==0) = 1;
    H_max = -sum(p_x1.*log2(p_x1));
end