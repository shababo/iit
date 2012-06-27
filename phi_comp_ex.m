function [phi prob prob_prod_MIP MIP] = phi_comp_ex(options,M,x0,x0_s,p,b_table,M_p,BRs,FRs)


% M_p: power set of M
% op_disp = 1;
% op_single = 1;


op_context = options(6);
op_empty = options(8);
op_min = options(9);
op_console = options(10);
op_volume = options(11);

N = length(M);
M_p{2^N} = []; % add empty set

if op_min == 0 % take sum of forward and backward
    phi_MIP = zeros(2^N-1,2^N-1);
    prob_cand = cell(2^N-1,2^N-1);
    prob_prod_MIP_cand = cell(2^N-1,2^N-1);
    MIP_cand = cell(2^N-1,2^N-1);
    
    if op_empty == 0
        i_max = 2^N-1;
    else
        i_max = 2^N;
    end
    
    for i = 1:i_max
        xp = M_p{i};
        for j=1: i_max
            xf = M_p{j};
            N_p = length(xp);
            N_f = length(xf);
            if N_p ~= 0 || N_f ~= 0
                if op_context == 0 % conservative
                    [phi_MIP(i,j) prob_cand{i,j} prob_prod_MIP_cand{i,j} MIP_cand{i,j}] ...
                        = phi_comp_bf(options,M,x0,xp,xf,x0_s,p,b_table,BRs,FRs);
                else % progressive
                    [phi_MIP(i,j) prob_cand{i,j} prob_prod_MIP_cand{i,j} MIP_cand{i,j}] ...
                        = phi_comp_bf(options,M,x0,xp,xf,x0_s,p,b_table);
                end
            end
        end
    end
    % exclusion principle
    [phi i j] = max2(phi_MIP,M_p);
    xp = M_p{i}; % perspective of the past
    xf = M_p{j}; % perpective of the future
    MIP = MIP_cand{i,j};
    prob = prob_cand{i,j};
    prob_prod_MIP = prob_prod_MIP_cand{i,j};
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE CURRENT SETTINGS TAKE US HERE    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
else % take minimum of forward and backward
    
    phi_MIP = zeros(2^N-1,2);
    prob_cand = cell(2^N-1,1);
    prob_prod_MIP_cand = cell(2^N-1,1);
    MIP_cand = cell(2^N-1,1);
    
    % for all purview demoniators (we use same set for forward and
    % backward) - can we pick these better, Larissa says yes!
    for i=1: 2^N-1
        x = M_p{i};
        [phi_MIP(i,:) prob_cand{i} prob_prod_MIP_cand{i} MIP_cand{i}] ...
            = phi_comp_bf(options,M,x0,x,x,x0_s,p,b_table,BRs,FRs);
    end
    
    % exlusion principle
    max_phi_MIP_bf = zeros(2,1); % backward and forward phi
    MIP = cell(2,2,2);
    prob = cell(2,1);
    prob_prod_MIP = cell(2,1);
    for bf = 1:2
        [max_phi_MIP_bf(bf) j_max] = max_ex(phi_MIP(:,bf),M_p);
        MIP(:,:,bf) = MIP_cand{j_max}(:,:,bf);
        prob{bf} = prob_cand{j_max}{bf};
        prob_prod_MIP{bf} = prob_prod_MIP_cand{j_max}{bf};
        if bf == 1
            xp = M_p{j_max};
        else
            xf = M_p{j_max};
        end
    end
    if (op_volume == 0)
        phi = min(max_phi_MIP_bf(1),max_phi_MIP_bf(2));
    else
        phi = max_phi_MIP_bf(1);
    end
end

%% imposing maxent on units outside of perspectives
if op_context == 0
    for i = 1:2
        if i == 1
            x = xp;
        else
            x = xf;
        end
        if length(x) ~= N
            prob{i} = expand_prob(prob{i},M,x);
            prob_prod_MIP{i} = expand_prob(prob_prod_MIP{i},M,x);
        end
    end
end

if op_console
    fprintf('Core concept: x0=%s xp=%s  xf=%s\n',mod_mat2str(x0),mod_mat2str(xp),mod_mat2str(xf));
    fprintf('phi=%f\n',phi);
end
% figure(1)
% subplot(1,2,1),imagesc(prob)
% subplot(1,2,2),imagesc(prob_prod_MIP)
% phi_MIP
% [phi i j] = max2(phi_MIP,M_p)
% pause;

end

function [X_max i_max j_max] = max2(X,M_p)
% exclusion principle: if the value is the same, take the bigger one
X_max = -Inf;
i_max = 1;
j_max = 1;
s_max = 0;
for i=1: size(X,1)
    for j=1: size(X,2)
        s = length(M_p{i}) + length(M_p{j});
        cond1 = X(i,j) > X_max;
        cond2 = X(i,j) == X_max && s>= s_max;
        if cond1 || cond2
            X_max = X(i,j);
            i_max = i;
            j_max = j;
            s_max = s;
        end
    end
end

end


function [X_max i_max] = max_ex(X,M_p)
% exclusion principle: if the value is the same, take the bigger one
X_max = -Inf;
i_max = 1;
s_max = 0;
for i=1: size(X,1)
    s = length(M_p{i});
    cond1 = X(i) > X_max;
    cond2 = X(i) == X_max && s>= s_max;
    if cond1 || cond2
        X_max = X(i);
        i_max = i;
        s_max = s;
    end
end

end