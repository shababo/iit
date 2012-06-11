function [phi MIP0 MIP1 prob prob1 prob2] = phi_comp(options,M,C_j,x1,p,b_table,x0_b1_o,x0_b2_o,N_b)

op_fb = options(1);
op_phi = options(2);
op_disp = options(3);
op_single = options(4);
op_ex = options(5);
% op_disp = 1;
% op_single = 1;

N = length(M);

% computing entropy in the whole system given C_j
if op_fb == 1
    % backward perspective
    [H prob] = partial_ent(op_fb,0,M,C_j,x1,p,b_table,M);
elseif op_fb == 0
    % forward perspective
    [H prob] = partial_ent(2,0,C_j,M,x1,p,b_table,M);
elseif op_fb == 3
    % forward and backward simultaneous perspective
    [H prob] = partial_ent(4,0,C_j,M,x1,p,b_table,M);
    % fprintf('%s\n',mat2str(prob))
end

i = length(C_j);

cond1 = i==1 && op_single == 0;
%% single perspective
if cond1 == 1 || length(M) == 1
     if op_fb == 0 || op_fb == 3
        [H_max p_x1] = forward_H(M,p,b_table,op_fb);
        phi = -H;
        for j=1: length(p_x1)
            if prob(j) ~=0 && p_x1(j) ~= 0
                phi = phi - prob(j)*log2(p_x1(j));
            elseif p_x1(j) == 0 && prob(j) ~= 0
                phi = Inf;
            end
        end
     elseif op_fb == 1
        H_max = N;
     end
     phi = H_max - H;
    
    if op_disp == 1 || op_disp == 2
        fprintf('phi=%f H_max=%f H=%f\n',phi,H_max,H);
    end
    MIP0 = M;
    MIP1 = C_j;
    prob1 = prob;
    prob2 = zeros(N,1);
else

%% more than one
x0_b1 = x0_b1_o;
x0_b2 = x0_b2_o;
if op_fb == 0 || op_fb == 3
    [x1_b1 x1_b2 N1_b] = bipartition(M,N); % partition of M
elseif op_fb == 1
    [x1_b1 x1_b2 N1_b] = bipartition(C_j,length(C_j)); % partition of C_j
end

phi_cand = zeros(N_b*N1_b,2);
H_vec = zeros(N_b*N1_b,2);
prob_vec = cell(N_b*N1_b,1);

for i_H = 1: N_b*N1_b
    j_b = rem(i_H,N1_b);
    if j_b == 0
        j_b = N1_b;
    end
    i_b = (i_H-j_b)/N1_b+1;
    x0_1 = x0_b1{i_b};
    x0_2 = x0_b2{i_b};
    x1_1 = x1_b1{j_b};
    x1_2 = x1_b2{j_b};
    
    if op_fb == 0 || op_fb == 3
        Norm = Normalization(x0_1,x1_1,x1_2,C_j);
    elseif op_fb == 1
        Norm = Normalization(x1_1,x0_1,x0_2,C_j);
    end  
    
    if Norm ~= 0
        prob_M = prob;
        [H1 prob1] = partial_ent(op_fb,op_phi,x0_1,x1_1,x1,p,b_table,M,C_j,prob_M);
        [H2 prob2] = partial_ent(op_fb,op_phi,x0_2,x1_2,x1,p,b_table,M,C_j,prob_M);
    else
        H1 = Inf;
        H2 = Inf;
        prob1 = [];
        prob2 = [];
    end
    
    
    % fprintf('%s %s %d\n',mat2str(x0_1),mat2str(x1_1),length(prob1));
    % fprintf('%s %s %d\n',mat2str(x0_2),mat2str(x1_2),length(prob2));    
    
    prob_vec{i_H} = [prob1; prob2];
    H_vec(i_H,1) = H1;
    H_vec(i_H,2) = H2;
    phi_cand(i_H,1) = (H1 + H2) - H;
    
    % Norm = min(length(x0_1),length(x0_2));
    % Norm = min(length(x0_1),length(x1_2)) + min(length(x0_2),length(x1_1));
    % Norm = 1; % no normalization
    phi_cand(i_H,2) = phi_cand(i_H,1)/Norm;
end

[min_norm_phi i_phi_min] = min(phi_cand(:,2));

j_b = rem(i_phi_min,N1_b);
if j_b == 0
    j_b = N1_b;
end
i_b = (i_phi_min-j_b)/N1_b+1;

phi = phi_cand(i_phi_min,1);
prob_MIP = prob_vec{i_phi_min};

MIP0 = x0_b1{i_b};
MIP1 = x1_b1{j_b};

if op_disp == 1 || op_disp == 2
    fprintf('phi=%f H=%f minH_sum=%f H1=%f H2=%f\n',phi,H,phi+H,H_vec(i_phi_min,1),H_vec(i_phi_min,2));
end

if op_fb == 0
    p1_length = length(MIP1);
elseif op_fb == 1
    p1_length = length(MIP0);
elseif op_fb == 3
    p1_length = 2*length(MIP1);
end

prob1 = prob_MIP(1:2^p1_length);
if 2^p1_length == length(prob_MIP)
     prob2 = [];
else
     prob2 = prob_MIP(2^p1_length+1:end);
end

% if length(prob1) ~= 0
% check_length = log2(length(prob1));
% end
% if length(prob2) ~= 0
%     check_length = check_length + log2(length(prob2));
% end
% if check_length ~= 2*N
%     fprintf('alert\n');
%     check_length
%     display(MIP0)
%     display(MIP1)
%     display(x0_b2{i_b})
%     display(x1_b2{j_b})
%     display(prob1')
%     display(prob2')
%     pause;
% end

end

end

function Norm = Normalization(s_1,t_1,t_2,source)
N_s = length(source);
N_s1 = 0;
for i=1: N_s
    N_s1 = N_s1 + sum(s_1 == source(i));
end
N_s2 = N_s - N_s1;
Norm = min(N_s1,length(t_2)) + min(N_s2,length(t_1));
end