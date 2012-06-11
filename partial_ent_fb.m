function [H prob_r] = partial_ent_fb(op_phi,xp_t,xp_out,x0_so,x0_in,x0_out,xf_t,source,p,b_table,prob_M)

% compute the conditional entropy of the present state xp_t and the future
% state xf_t given the current state x0_so, H(xp_t,xf_t|x0_so)
% source: data of x0
% p: probability matrix in the whole system

T = 1;
N = size(p,2);

%% simultaneous computation of backward and forward
if isempty(x1_b) == 1
    H = 0;
    prob_r = [];
else
    %% perspective
    if op_phi == 1 && nargin < 10
        x0_so = C_j; % target
        x0_in = pick_rest(M,C_j);
        x0_out = pick_rest(1:N,M);
        x1_in = M; % source
        x1_out = pick_rest(1:N,M);
        prob_M =   partial_prob_comp_bf(T,xp_t,xp_out,x0_so,x0_in,x0_out,xf_t,source,p,b_table);
    end
    
    if op_fb == 3
        %% partition
        x0_so = [];
        for j=1: length(C_j)
            x0_so =  [x0_so x0_b(x0_b==C_j(j))];
        end
        x0_in = pick_rest(x0_b,x0_so);
        x0_out = pick_rest(1:N,x0_b);
    elseif op_fb == 4
        %% perspective
        x0_so = x0_b;
        x0_in = pick_rest(M,x0_so);
        x0_out = pick_rest(1:N,M);
    end
    x1_in = x1_b;
    x1_out = pick_rest(1:N,x1_b);
    % fprintf('x0_so=%s x0_in=%s x0_out=%s\n',mat2str(x0_so),mat2str(x0_in),mat2str(x0_out))
    prob = patial_prob_comp_bf(T,x0_so,x0_in,x0_out,x1_in,x1_out,source,p,b_table);
    
    % subplot(1,3,1),imagesc(prob);
    % subplot(1,3,2),imagesc(prob_M);
    prob = prob(:);
    
    if op_phi == 1
        prob2 = p_marginalization(M,x1_b,b_table,prob_M,op_fb);
        % subplot(1,3,3),imagesc(prob2);
        prob2 = prob2(:);
    end
    
    prob_r = prob;
    
    if sum(prob) ~= 0
        prob(prob==0) = 1; % avoid log0 when computing entropy
        if op_phi == 0
            H = -sum(prob.*log2(prob)); % conditional entropy H(x0_b|x1_b)
        else
            H = -sum(prob2.*log2(prob));
        end
    else
        H = 0;
    end
    
end