clear all;
close all;

tic;

%% options

%put options in array... TODO: this should be a struct and maybe also
%global

op_complex = 1;  % 0: only consider the whole system %1: find the complex
op_parallel = 0; % 1: parallel computing, 0: not
op_network = 4; % 1: random with self-connectivity, 2: random without self-connectivity, 3: modular, 4: logic gates, 0: some given connectivity matrix
op_TPM = 0; % 1: load TPM
op_ave = 0; % 0: use a specific current state 1: average over all possible current states
op_disp = 0; % 0: No figures, 1: only complex, 2: complex and whole system, 3: all figures
op_context = 0; % 0: conservative 1: progressive
op_empty = 1; % 0: excluding empty set in the past and the future 1: including empty set 
op_min = 1; % conservative only 0: phi is the sum of phi backward and phi forward (simulataneous partition)
                     % 1: phi is the minimum of phi_b and phi_f (separate partition)
                     
%% inactive options, which are not used anymore
op_fb = 3; % 0: forward repertoire, 1: backward repertoire 2: both separately 3: both simultaneously
op_phi = 1; % two versions of small phi 0:Difference of entropy, 1:KL-divergence
op_whole = 0; % KLD is computed in 0: small system 1: whole system (previous version)


options = [op_fb op_phi op_disp 1 1 op_context op_whole op_empty op_min];

save options options

%% define the connectivty of the network
N = 3; % Number of elements in the network
Na = 3; % Number of afferent connections

% current state
if op_ave == 0
    current_state = zeros(N,1); % all OFF
    % current_state = ones(N,1); % all ON
    % current_state = [0 1 0 1]';
    z_max = 1;
else
    z_max = 2^N;
end

J = zeros(N,N); % connectivity matrix
if op_network == 1     % random network with self-connectivity
    
    for i=1: N % for each element
        x = 1:N;
        y = randsample(N,Na); %Na samples from list N w/o replacement, if N=Na, then this is equiv to randperm(Na)
        z = x(y);
        J(i,z) = 1;
        % J(i,z) = 0.5/Na;
    end
    
elseif op_network == 2
    % random network without self-connectivity
    for i=1: N
        x = 1: N;
        x(i) = [];
        y = randsample(N-1,Na);
        z = x(y);
        J(i,z) = 1;
        % J(i,z) = 0.5/Na;
    end
elseif op_network == 3
    % modular network
    for i=1: N/2
        i1 = 2*i-1;
        i2 = 2*i;
        J(i1:i2,i1:i2) = (N-1)*[0 1; 1 0];
    end
elseif op_network == 4
    %% logic gates
    % logic type: 1-> AND, 2-> OR, 3-> XOR, 4 -> COPY, 5-> NOT, 6 -> NULL
    logic_type = zeros(N,1);
    % 2 COPY
%     logic_type(1) = 4;
%     logic_type(2) = 4;
%     J(1,2) = 1;
%     J(2,1) = 1;
    
    % 1 XOR, 1 AND, 2 ORs
%     logic_type(1) = 3;
%     logic_type(2) = 1;
%     logic_type(3) = 2;
%     logic_type(4) = 2;
%     J(1,[3 4]) = 1;
%     J(2,[3 4]) = 1;
%     J(3,[1 4]) = 1;
%     J(4,[2 3]) = 1;

    % 1 XOR, 2 OR
    logic_type(1) = 3;
    logic_type(2) = 2;
    logic_type(3) = 2;
    J(1,[2 3]) = 1;
    J(2,[1 3]) = 1;
    J(3,[1 2]) = 1;

% 1 AND, 2 NULL
%     logic_type(1) = 1;
%     logic_type(2) = 6;
%     logic_type(3) = 6;
%     J(1,[2 3]) = 1;
%     J(2, []) = 1;
%     J(3, []) = 1;

    % 1 AND, 2 COPY
%     logic_type(1) = 4;
%     logic_type(2) = 1;
%     logic_type(3) = 4;
%     J(1,[2]) = 1;
%     J(2, [1 3]) = 1;
%     J(3, [2]) = 1;

% 1 AND, 1 COPY, 2 NULL
%     logic_type(1) = 1;
%     logic_type(2) = 6;
%     logic_type(3) = 6;
%     logic_type(4) = 4;
%     J(1,[2 3]) = 1;
%     % J(2,[3 4]) = 1;
%     % J(3,[1 4]) = 1;
%     J(4,[1]) = 1;
else
    load J_fix;
end

%% compute the transition probability matrix
p_x0 = zeros(2^N,N);
if op_network == 4
    % logic gates
    for k=1: 2^N
        x0 = trans2(k-1,N);
        for i=1: N
            i_vec = find(J(i,:)==1);
            input_vec = x0(i_vec);
            p_x0(k,i) = logic_gates(input_vec,logic_type(i));
        end
    end
else
    % sigmoid function
    % a: threshold T: the level of noise
    a = floor(Na/2) + 0.5; % majority vote
    % a = N;
    % a = max(max(J)) - 0.1; % generalized OR
    T = 0.01; % T=0: no noise
    for i=1: 2^N
        x0 = trans2(i-1,N);
        p_x0(i,:) = (1+tanh((J*x0-a)/T))/2; % probability of turning on given the data x0
        % p_x0(i,:) = J*x0;
    end
end

%% p(t=1|t=0)

% This section computes the forward probability of each state given the
% current state

p1 = ones(2^N,2^N); % p(x0,x1) = p(x1|x0)

for i=1: 2^N
    for j=1: 2^N
        x1 = trans2(j-1,N);
        for k=1: N
            if x1(k) == 1
                p1(i,j) = p1(i,j)*p_x0(i,k);
            else
                p1(i,j) = p1(i,j)*(1-p_x0(i,k));
            end
        end
    end
end

% p = p1; % practical version
if op_TPM == 0
    p = p_x0; % original version
else
    load TPM;
    N = size(p,2);
end

% if op_disp > 1
%     figure(100)
%     subplot(1,2,1),imagesc(J)
%     colormap('gray')
%     subplot(1,2,2),imagesc(p)
% end
% pause;

f = zeros(N,1); % averaged firing rates
for i=1: N
    f(i) = sum(p_x0(:,i))/2^N;
end
avef = sum(f)/N;
fprintf('avef=%f\n',avef);

p_x1 = zeros(2^N,1);
for i=1: 2^N
    x1 = trans2(i-1,N);
    for j=1: 2^N
        p_temp = 1;
        for k=1: N
            if x1(k) == 1
                p_temp = p_temp*p_x0(j,k);
            else
                p_temp = p_temp*(1-p_x0(j,k));
            end
        end
        p_x1(i) = p_x1(i) + p_temp;
    end
end
p_x1 = p_x1/2^N;

%% binary table
b_table = cell(2^N,N);
states = zeros(N,2^N);
for i=1: N
    for j=1: 2^i
        b_table{j,i} = trans2(j-1,i);
        if i== N
            states(:,j) = trans2(j-1,i);
        end
    end
end

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
%     H_max = N;
% end

Big_phi_st = zeros(2^N,1);
Big_phi_MIP_st = zeros(2^N,1);

for z=1: z_max
    if op_ave == 0
        x1 = current_state;
    else
        x1 = trans2(z-1,N);
    end
    fprintf('x1=%s\n',mat2str(x1));
    
    % partial_prob_comp(partition, partition, state, prob_matrix, binary
    % table, op_fb
    check_prob = partial_prob_comp(1:N,1:N,x1,p,b_table,1);
    state_check = sum(check_prob);
    if state_check == 0
        fprintf('This state cannot be realized!\n')
        Big_phi_st(z) = 0;
        Big_phi_MIP_st(z) = 0;
    else
        if op_complex == 0 % only consider whole system
            if op_fb == 2
                options(1) = 0; [Big_phi_f phi_f prob_cell_f] = big_phi_comp(1:N,x1,p,b_table,options);
                options(1) = 1; [Big_phi_b phi_b prob_cell_b] = big_phi_comp(1:N,x1,p,b_table,options);
                Big_phi = Big_phi_f + Big_phi_b;
            elseif op_fb == 3 % THIS IS THE ONLY ONE WE DO NOW? BOTH FORWARD AND BACKWARD SIMULTANEOUSLY
                M = 1:N;
                if op_context == 0 % THIS IS REDUNDANT - WE COULDN'T BE HERE IF THIS WEREN'T TRUE!!
                    [BRs FRs] = comp_pers(x1,p,b_table,options);
                    [Big_phi phi prob_cell MIP prob_cell2] = big_phi_comp_fb(M,x1,p,b_table,options,BRs,FRs);
                else
                    [Big_phi phi prob_cell MIP prob_cell2] = big_phi_comp_fb(M,x1,p,b_table,options);
                end
            else
                [Big_phi phi prob_cell] = big_phi_comp(1:N,x1,p,b_table,options);
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
        Big_phi_ave = sum(Big_phi_st)/2^N;
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