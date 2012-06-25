% THIS SCRIPT COMPUTES BIG PHI OVER A BINARY NETWORK THAT IS CURRENTLY SET
% TO BE A SPECIFIC LOGICAL NETWORK

clear all;
close all;

tic;

%% options

%put options in array... TODO: this should be a struct and maybe also
%global

global options

options.complex = 1; % 0: only consider the whole system %1: find the complex
options.parallel = 0; % 1: parallel computing, 0: not
options.network_type = 4; % 1: random with self-connectivity, 2: random without self-connectivity, 3: modular, 4: logic gates, 0: some given connectivity matrix
options.tpm = 0;  % 1: load TPM
options.states = 1;  % 0: use a specific current state 1: average over all possible current states
options.display = 0; % 0: No figures, 1: only complex, 2: complex and whole system, 3: all figures
options.context = 0;  % 0: conservative 1: progressive <-- probably lose progressive
options.minimum = 1; % conservative only 0: phi is the sum of phi backward and phi forward (simulataneous partition)
                     % 1: phi is the minimum of phi_b and phi_f (separate
                     % partition) <--- we'll only phi in min
                     
%% inactive options, which are not used anymore
op_fb = 3; % 0: forward repertoire, 1: backward repertoire 2: both separately 3: both simultaneously
op_phi = 1; % two versions of small phi 0:Difference of entropy, 1:KL-divergence
op_whole = 0; % KLD is computed in 0: small system 1: whole system (previous version)



save options options

%% define the connectivty of the network
num_elements = 3; % Number of elements in the network
num_afferents = 3; % Number of afferent connections
num_states = 2^num_elements; % number of total states in the network


% current state
if options.states == 0 % we only check a single state
    current_state = zeros(num_elements,1); % all OFF
    % current_state = ones(num_elements,1); % all ON
    % current_state = [0 1 0 1]';
    num_states_of_interest = 1;
else
    num_states_of_interest = num_states;
end

connectivity_matrix = zeros(num_elements,num_elements); % connectivity matrix

if options.network_type == 1     % random network with self-connectivity
    
    for i=1: num_elements % for each element
        x = 1:num_elements;
        y = randsample(num_elements,num_afferents); %num_afferents samples from list num_elements w/o replacement, if num_elements=num_afferents, then this is equiv to randperm(num_afferents)
        z = x(y);
        connectivity_matrix(i,z) = 1;
        % connectivity_matrix(i,z) = 0.5/num_afferents;
    end
    
elseif options.network_type == 2
    % random network without self-connectivity
    for i=1: num_elements
        x = 1: num_elements;
        x(i) = [];
        y = randsample(num_elements-1,num_afferents);
        z = x(y);
        connectivity_matrix(i,z) = 1;
        % connectivity_matrix(i,z) = 0.5/num_afferents;
    end
elseif options.network_type == 3
    % modular network
    for i=1: num_elements/2
        i1 = 2*i-1;
        i2 = 2*i;
        connectivity_matrix(i1:i2,i1:i2) = (num_elements-1)*[0 1; 1 0];
    end
elseif options.network_type == 4
    %% logic gates
    % logic type: 1-> AND, 2-> OR, 3-> XOR, 4 -> COPY, 5-> NOT, 6 -> NULL
    logic_type = zeros(num_elements,1);
    % 2 COPY
%     logic_type(1) = 4;
%     logic_type(2) = 4;
%     connectivity_matrix(1,2) = 1;
%     connectivity_matrix(2,1) = 1;
    
    % 1 XOR, 1 AND, 2 ORs
%     logic_type(1) = 3;
%     logic_type(2) = 1;
%     logic_type(3) = 2;
%     logic_type(4) = 2;
%     connectivity_matrix(1,[3 4]) = 1;
%     connectivity_matrix(2,[3 4]) = 1;
%     connectivity_matrix(3,[1 4]) = 1;
%     connectivity_matrix(4,[2 3]) = 1;

    % 1 XOR, 2 OR
    logic_type(1) = 3;
    logic_type(2) = 2;
    logic_type(3) = 2;
    connectivity_matrix(1,[2 3]) = 1;
    connectivity_matrix(2,[1 3]) = 1;
    connectivity_matrix(3,[1 2]) = 1;

% 1 AND, 2 NULL
%     logic_type(1) = 1;
%     logic_type(2) = 6;
%     logic_type(3) = 6;
%     connectivity_matrix(1,[2 3]) = 1;
%     connectivity_matrix(2, []) = 1;
%     connectivity_matrix(3, []) = 1;

    % 1 AND, 2 COPY
%     logic_type(1) = 4;
%     logic_type(2) = 1;
%     logic_type(3) = 4;
%     connectivity_matrix(1,[2]) = 1;
%     connectivity_matrix(2, [1 3]) = 1;
%     connectivity_matrix(3, [2]) = 1;

% 1 AND, 1 COPY, 2 NULL
%     logic_type(1) = 1;
%     logic_type(2) = 6;
%     logic_type(3) = 6;
%     logic_type(4) = 4;
%     connectivity_matrix(1,[2 3]) = 1;
%     % connectivity_matrix(2,[3 4]) = 1;
%     % connectivity_matrix(3,[1 4]) = 1;
%     connectivity_matrix(4,[1]) = 1;
else
    load J_fix;
end

%% binary table
b_table = cell(num_states,num_elements);

for i= 1:num_elements
    for j= 1:2^i
        b_table{j,i} = trans2(j-1,i);
    end
end

%% compute the transition probability matrix

p_x0 = zeros(num_states,num_elements);

if options.network_type == 4    % logic gates
    % for each state
    for k = 1:num_states
        x0 = trans2(k-1,num_elements);
        for i = 1:num_elements
            i_vec = find(connectivity_matrix(i,:)==1);
            input_vec = x0(i_vec);
            p_x0(k,i) = logic_gates(input_vec,logic_type(i));
        end
    end
else
    % sigmoid function
    % a: threshold T: the level of noise
    a = floor(num_afferents/2) + 0.5; % majority vote
    % a = num_elements;
    % a = max(max(connectivity_matrix)) - 0.1; % generalized OR
    T = 0.01; % T=0: no noise
    for i=1: num_states
        x0 = trans2(i-1,num_elements);
        p_x0(i,:) = (1+tanh((connectivity_matrix*x0-a)/T))/2; % probability of turning on given the data x0
        % p_x0(i,:) = connectivity_matrix*x0;
    end
end



%% p(t=1|t=0)

% This section computes the forward probability of each state given the
% current state

p1 = ones(num_states,num_states); % p(x0,x1) = p(x1|x0)

for i=1: num_states
    for connectivity_matrix=1: num_states
        x1 = trans2(connectivity_matrix-1,num_elements);
        for k=1: num_elements
            if x1(k) == 1
                p1(i,connectivity_matrix) = p1(i,connectivity_matrix)*p_x0(i,k);
            else
                p1(i,connectivity_matrix) = p1(i,connectivity_matrix)*(1-p_x0(i,k));
            end
        end
    end
end

% p = p1; % practical version
if op_TPM == 0
    p = p_x0; % original version
else
    load TPM;
    num_elements = size(p,2);
end

% if op_disp > 1
%     figure(100)
%     subplot(1,2,1),imagesc(connectivity_matrix)
%     colormap('gray')
%     subplot(1,2,2),imagesc(p)
% end
% pause;

f = zeros(num_elements,1); % averaged firing rates - THIS ASSUMES THAT ALL STATES ARE EQUALLY REPRESENTED DURING NORMAL FIRING PATTERNS!!
for i = 1:num_elements
    f(i) = sum(p_x0(:,i))/num_states;
end
avef = sum(f)/num_elements;

fprintf('avef=%f\n',avef);

p_x1 = zeros(num_states,1);

%%%%%%% THE HUGE FOR-LOOPS BELOW CAN ACTUALLY BE COMPUTED BY SUMMING OVER
%%%%%%% THE COLUMNS OF p1 AND DIVIDING BY THE # OF COLS... ALSO, p1 IS
%%%%%%% NEVER USED IN THIS CODE

for i=1: num_states % current state of interest
    x1 = trans2(i-1,num_elements);
    for j = 1:num_states % previous state
        p_temp = 1;
        for k=1: num_elements % checks over each element
            % p_temp == 1 IFF x1 == p_x0(j,:)
            % basically we are counting how many times we go into state x1
            if x1(k) == 1 
                p_temp = p_temp*p_x0(j,k);
            else
                p_temp = p_temp*(1-p_x0(j,k));
            end
        end
        p_x1(i) = p_x1(i) + p_temp;
    end
end
p_x1 = p_x1/num_states; % normalize to probability



%% parallel computing
isOpen = matlabpool('size');
if  isOpen == 0 && op_parallel > 0
%     s = ['matlabpool ' int2str(op_parallel)];
%     eval(s);
    matlabpool;
end

%% compute the big-phi
% if op_fb == 0
%     p = p'; % forward repertoire
% else
%     H_max = num_elements;
% end

Big_phi_st = zeros(num_states,1);
Big_phi_MIP_st = zeros(num_states,1);

% for each "current state"
for z = 1:num_states_of_interest
    if op_ave == 0
        x1 = current_state;
    else
        x1 = trans2(z-1,num_elements);
    end
    fprintf('x1=%s\n',mat2str(x1));
    
    % partial_prob_comp(partition, partition, state, prob_matrix, binary
    % table, op_fb
    check_prob = partial_prob_comp(1:num_elements,1:num_elements,x1,p,b_table,1);
    state_check = sum(check_prob);
    if state_check == 0
        fprintf('This state cannot be realized!\n')
        Big_phi_st(z) = 0;
        Big_phi_MIP_st(z) = 0;
    else
        if op_complex == 0 % only consider whole system
            if op_fb == 2
                options(1) = 0; [Big_phi_f phi_f prob_cell_f] = big_phi_comp(1:num_elements,x1,p,b_table,options);
                options(1) = 1; [Big_phi_b phi_b prob_cell_b] = big_phi_comp(1:num_elements,x1,p,b_table,options);
                Big_phi = Big_phi_f + Big_phi_b;
            elseif op_fb == 3 % THIS IS THE ONLY ONE WE DO NOW? BOTH FORWARD AND BACKWARD SIMULTANEOUSLY
                M = 1:num_elements;
                if op_context == 0
                    [BRs FRs] = comp_pers(x1,p,b_table,options);
                    [Big_phi phi prob_cell MIP prob_cell2] = big_phi_comp_fb(M,x1,p,b_table,options,BRs,FRs);
                else
                    [Big_phi phi prob_cell MIP prob_cell2] = big_phi_comp_fb(M,x1,p,b_table,options);
                end
            else
                [Big_phi phi prob_cell] = big_phi_comp(1:num_elements,x1,p,b_table,options);
            end
            Big_phi_st(z) = Big_phi;
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % THE CURRENT SETTINGS TAKE US HERE    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else % find the complex
            [Big_phi_MIP MIP Big_phi_M IRR_phi IRR_REP IRR_MIP M_IRR prob_M phi_M MIP_M] ...
                = big_phi_complex(x1,p,b_table,options);
            
            if op_fb == 2
                % subindex b means backward and f means forward
                IRR_phi_b = IRR_phi{1};
                IRR_phi_f = IRR_phi{2};
                IRR_REP_b = IRR_REP{1};
                IRR_REP_f = IRR_REP{2};
                M_IRR_b = M_IRR{1};
                M_IRR_f = M_IRR{2};
                % state dependent big phi and big phi MIP
                Big_phi_st(z) = sum(IRR_phi_b)+sum(IRR_phi_f);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % THE CURRENT SETTINGS TAKE US HERE    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                Big_phi_st(z) = sum(IRR_phi);
            end
            Big_phi_MIP_st(z) = Big_phi_MIP;
        end
    end
    
    % pause;
end

if op_ave == 1
    if op_fb == 0
        Big_phi_ave = sum(Big_phi_st)/num_states;
    else
        Big_phi_ave = sum(p_x1 .* Big_phi_st); %hmmm... interesting, weighted ave
    end
    fprintf('Big_phi_ave=%f\n',Big_phi_ave);
end

op_close = 0;
isOpen = matlabpool('size');
if isOpen > 0 && op_close == 1
    matlabpool close;
end

toc;