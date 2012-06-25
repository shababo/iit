function [Big_phi phi prob_cell MIP prob_cell2] = big_phi_comp_fb(M,x0_s,p,b_table,options,BRs,FRs)

%%  compute big phi for a subset, M
% M: a subset of the whole system (can be the whole system itself)
% x0_s: given data about the current state
% p: transition probability matrix in the whole system (p(x1|x0))

% THE FINAL TWO ARGS ARE OPTIONS, IF THEY ARE THERE THEN WE ARE DOING
% CONSERVATIVE, OTHERWISE PROGRESSIVE...

N = length(M);

op_fb = options(1);
op_phi = options(2);
op_disp = options(3);
op_single = options(4);
op_context = options(6);
op_min = options(9);

op_whole = 1;

if op_disp == 0 || op_disp == 1
    disp_flag = 0;
elseif op_disp == 2 && N ~= log2(size(p,1)) % or we are not dealing with the whole complex
    disp_flag = 0;
else
    disp_flag = 1;
    if op_fb == 0
        fig_p = 200;
    else
        fig_p = 100;
    end
end

%% x0 data


C_x0 = cell(2^N-1,1);
k = 1;
for i = 1:N % can this be done in one for-loop over k = 1:2^N-1 ?
    C = nchoosek(M,i); % create a matrix of combinations of M of size i
    N_C = size(C,1);
    % fprintf('i=%d N_c=%d\n',i,N_C);
    for j = 1:N_C % for all combos of size i
        x0 = C(j,:); % pick a combination
        C_x0{k} = x0;% store combo
        k = k + 1;
    end
end
M_p = C_x0; % power set of M

MIP = cell(2^N-1,1); % MIP in the past, the present, and the future
phi = zeros(2^N-1,1); % small phis

prob = cell(2^N-1,1); % transition repertoire
prob_prod = cell(2^N-1,1); % partitioned transition repertoire

%% computing small phis
parfor ci=1: 2^N-1
    x0 = C_x0{ci}; % given data of x0
    if op_disp ~= 2
        fprintf('C=%s\n',mod_mat2str(x0));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THE CURRENT SETTINGS TAKE US HERE    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if op_context == 0
        [phi(ci) prob{ci} prob_prod{ci} MIP{ci}] ...
        =  phi_comp_ex(options,M,x0,x0_s,p,b_table, M_p, BRs, FRs);
    else
        [phi(ci) prob{ci} prob_prod{ci} MIP{ci}] ...
            =  phi_comp_ex(options,M,x0,x0_s,p,b_table, M_p);
    end
end

prob_cell = cell(2,1);
prob_cell{1} = prob;
prob_cell{2} = prob_prod;

phi_m = zeros(N,3); % cumulative sum

% PRETTY SURE THIS CAN JUST BE DONE WITH A SUM() CALL
for i_C=1: 2^N-1
    C = C_x0{i_C};
    i = length(C);
    phi_m(i,1) = phi_m(i,1) + phi(i_C);
    phi_m(i,2) = phi_m(i,2) + phi(i_C)/nchoosek(N,i);
end

% THIS SEEMS LIKE IT CAN BE DONE A SMARTER WAY
for i=1: N
    if i > 1
        phi_m(i,3) = phi_m(i-1,3) + phi_m(i,1);
    else
        phi_m(i,3) = phi_m(i,1);
    end
end

Big_phi = phi_m(end,3);

%% display
for i_C=1: 2^N-1
    C = C_x0{i_C};
    i = length(C);
    
    if i > 1 || op_single == 1        
        [string_p string] = make_title_fb(MIP{i_C},op_context,op_min);
        fprintf('C=%s: Core Concept: %s\n',mod_mat2str(C),string{3});
        fprintf('Partition: %s',string_p{3});
    end
    
    if abs(phi(i_C)) < 10^-8
        phi(i_C) = 0;
    end
    fprintf(' phi%d=%f\n',i,phi(i_C));
    % prob{i_C}
    % prob_prod{i_C}
end

prob_cell2 = cell(2^N-1,2);
for i=1: 2^N-1
    prob_cell2{i,1} = prob{i};
    prob_cell2{i,2} = prob_prod{i};
end
if disp_flag == 1
    if size(p,2) == N
        fprintf('\n')
        fprintf('------------------------------------------------------------------------\n')
        fprintf('Whole system\n')
        fprintf('Core concepts: MIP: Small phi\n')
    end
    plot_REP(prob_cell2,phi,MIP,100, M, op_context, op_min);
end

for i=1: N
    fprintf('%d: phi_cum=%f phi_sum=%f phi_mean=%f\n',i,phi_m(i,3),phi_m(i,1),phi_m(i,2));
end

fprintf('Big phi=%f\n',Big_phi);

end
