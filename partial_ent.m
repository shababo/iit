function [H prob_r] = partial_ent(op_fb,op_phi,x0_b,x1_b,source,p,b_table,M,C_j,prob_M)

% compute the conditional entropy, H(x0_b|x1_b), given a partition of x0 and x1
% x1: data, x0_b: partition of x0, x1_b: partition of x1
% p: probability matrix in the whole system
% C_j: subset of x1

T = 1;
N = size(p,2);

%% simultaneous computation of backward and forward
if op_fb >= 3
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
    
    
%% separate computation of backward and forward
else

if op_fb == 0 && isempty(x1_b) == 1
    H = 0;
    prob_r = [];
elseif op_fb == 1 && isempty(x0_b) == 1
    H = 0;
    prob_r = [];
else

if exist('C_j','var') == 0
    C_j = [];
end

%% perspective
if op_phi == 1
    if op_fb == 1
        x0_so = M; % target
        x0_in = [];
        x0_out = pick_rest(1:N,x0_so);
        x1_tar = C_j; % source
        % prob_M = partial_prob_comp2(x0_tar,x0_in,x0_out,x1_so,source,p,b_table,op_fb);
        % prob_M = partial_prob_comp(M,C_j,x1,p,b_table,op_fb,M,C_j);
    else
        x0_so = C_j; % source
        x0_in = pick_rest(M,C_j);
        x0_out = pick_rest(1:N,M);
        x1_tar = M; % target
        % prob_M = partial_prob_comp2(x0_so,x0_in,x0_out,x1_tar,source,p,b_table,op_fb);
        % prob_M = partial_prob_comp(C_j,M,x1,p,b_table,2,M,C_j);
    end
    prob_M = partial_prob_comp_time(T,x0_so,x0_in,x0_out,x1_tar,source,p,b_table,op_fb);
end


%% partition (op_fb=0,1) or perspective (op_fb=2)
% prob = partial_prob_comp(x0_b,x1_b,x1,p,b_table,op_fb,M,C_j);
if op_fb == 1
    % perspective or partition (no differnece)
    x0_so = x0_b; % target
    x0_in = [];
    x0_out = pick_rest(1:N,x0_so);
    x1_tar = x1_b; % soruce
    % prob = partial_prob_comp2(x0_tar,x0_in,x0_out,x1_so,source,p,b_table,op_fb);
elseif op_fb == 0
    % partition
    x0_so = [];
    for j=1: length(C_j)
        x0_so =  [x0_so x0_b(x0_b==C_j(j))];
    end
    x0_in = pick_rest(x0_b,x0_so);
    x0_out = pick_rest(1:N,x0_b);    
    x1_tar = x1_b; % target
    % prob = partial_prob_comp2(x0_so,x0_in,x0_out,x1_tar,source,p,b_table,op_fb);
elseif op_fb == 2
    % perspective
    x0_so = x0_b;
    x0_in = pick_rest(M,x0_so);
    x0_out = pick_rest(1:N,M);
    x1_tar = x1_b; % target
    % prob =  partial_prob_comp2(x0_so,x0_in,x0_out,x1_tar,source,p,b_table,op_fb);   
end
prob = partial_prob_comp_time(T,x0_so,x0_in,x0_out,x1_tar,source,p,b_table,op_fb);

if op_phi == 1
    if op_fb == 1
        prob2 = p_marginalization(M,x0_b,b_table,prob_M);
    else
        prob2 = p_marginalization(M,x1_b,b_table,prob_M);
    end
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

end