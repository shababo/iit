clear all;
close all;

fprintf('\nRunning...\n\n');
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
op_console = 0; % 0: limited console output, 1: full console output
op_big_phi = 4; % 0 = big_phi is sum of small phi, 1 = big phi is volume best on EMD/best-packing, 2 = ave of H-difference b/w parent/child in hasse diagram
                     % 3 = look for pairwise distances less than 2*radiusas marker for overlap % 4 = take distance b/w whole and partitioned concepts and add to sum of change in
                     % small phi for each concept - FOR THIS OPTION PLEASE
                     % SET OP_BIG_PHI_DIST
op_big_phi_dist = 1; % 0 = KLD, 1 = EMD;
op_small_phi = 0; % 0 = use KLD between distributions, 1 = use EMD between distributions
op_sum = 0; % 0 = compute big_phi_mip based on expanding parts into the space of the whole, 1 = just take whole minus sum of parts
op_normalize_big_phi = 1; % 0 = don't normalize big_phi, 1 = normalize big_phi but choose non-norm value, 2 = normalize but choose norm value
op_normalize_small_phi = 1; % 0 = don't normalize small_phi, 1 = normalize small_phi but choose non-norm value, 2 = normalize but choose norm value

if (op_big_phi == 4)
    op_complex = 1;
end


%% inactive options, which are not used anymore
op_fb = 3; % 0: forward repertoire, 1: backward repertoire 2: both separately 3: both simultaneously
op_phi = 1; % two versions of small phi 0:Difference of entropy, 1:KL-divergence
op_whole = 0; % KLD is computed in 0: small system 1: whole system (previous version)

global grain, global noise;
grain = 100;
noise = 0.0; % 0 <= noise <= .5, this noise is applied to all nodes on their output essentially adding uncertaintity
global BRs, global FRs


fprintf('Noise Level = %f\n',noise);

options = [op_fb op_phi op_disp 1 1 op_context op_whole op_empty op_min op_console op_big_phi op_sum...
           op_normalize_big_phi op_normalize_small_phi op_complex op_small_phi op_big_phi_dist];

save options options

%% define the connectivty of the network
N = 3; % Number of elements in the network %!!!!!!!!!!!! CAN WE MAKE THIS DEPENDENT?
Na = 3; % Number of afferent connections

BRs = cell(2^N,2^N); % backward repertoire
FRs = cell(2^N,2^N); % forward repertoire


% current state
if op_ave == 0
%     current_state = zeros(N,1); % all OFF
    current_state = ones(N,1); % all ON
%     current_state = [1 1 0 0 0 0 0 0]';
%     current_state = [1 0]';
    z_max = 1;
else
    z_max = 2^N;
end

global J

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
    % logic type: 1-> AND, 2-> OR, 3-> XOR, 4 -> COPY, 5-> NOT, 6 -> NULL,
    % 7 -> MAJORITY, 8-> MINORITY, 9 -> PARITY
    logic_type = zeros(N,1);
    % 2 COPY
%     logic_type(1) = 4;
%     logic_type(2) = 4;
%     J(1,2) = 1;
%     J(2,1) = 1;

% THIS IS MY 5 NODE NETWORK!!!!

    % 1 XOR, 2 AND, 2 ORs
%     logic_type(1) = 3;
%     logic_type(2) = 1;
%     logic_type(3) = 2;
%     logic_type(4) = 2;
%     logic_type(5) = 1;
%     J(1,[3 4]) = 1;
%     J(2,[3 4]) = 1;
%     J(3,[1 4]) = 1;
%     J(4,[2 3]) = 1;
%     J(5,[2 4]) = 1;

% ---------------------------------------------------------------
% NORMALIZATION EXPERIMENTS

% 1) PARITY, MAJORITY, OR
    logic_type(1) = 2;
    logic_type(2) = 7;
    logic_type(3) = 9;
    J(1,[2 3]) = 1;
    J(2,[1 2 3]) = 1;
    J(3,[1 2 3]) = 1;
    
% 2) HOMOGENOUS AND
%     logic_type(1) = 1;
%     logic_type(2) = 1;
%     logic_type(3) = 1;
%     logic_type(4) = 1;
%     logic_type(5) = 1;
%     J(1,[1 2 3 4 5]) = 1;
%     J(2,[1 2 3 4 5]) = 1;
%     J(3,[1 2 3 4 5]) = 1;
%     J(4,[1 2 3 4 5]) = 1;
%     J(5,[1 2 3 4 5]) = 1;  
    
%     logic_type(1) = 1;
%     logic_type(2) = 1;
%     logic_type(3) = 1;
%     J(1,[1 2 3]) = 1;
%     J(2,[1 2 3]) = 1;
%     J(3,[1 2 3]) = 1;
    
% 3) COMPLEX NETWORK BASED ON PARITY/MAJORITY READING FOUR NODES BUT ADDED
% BACK CONNECTIONS SO WE GET ACTUAL VALUES FOR PHI/BIG_PHI/COMPLEX
%     logic_type(1) = 4;
%     logic_type(2) = 2;
%     logic_type(3) = 1;
%     logic_type(4) = 4;
%     logic_type(5) = 9;
%     logic_type(6) = 8;
%     J(1,5) = 1;
%     J(2,[5 6]) = 1;
%     J(3,[5 6]) = 1;
%     J(4,6) = 1;
%     J(5,[1 2 3 4]) = 1;
%     J(6,[1 2 3 4]) = 1;

% 4) MAJORITY - NON-HOMOGENOUS (EACH NODE HAS 3 AFFERENTS/EFFERENTS
%     logic_type(1) = 7;
%     logic_type(2) = 7;
%     logic_type(3) = 7;
%     logic_type(4) = 7;
%     logic_type(5) = 7;
%     J(1,[1 2 3]) = 1;
%     J(2,[3 4 5]) = 1;
%     J(3,[1 2 5]) = 1;
%     J(4,[1 4 5]) = 1;
%     J(5,[2 3 4]) = 1;

% 5) MODULAR COPIES (CLASSIC EXAMPLE OF FULL SET NOT BEING THE COMPLEX)
%     logic_type(1) = 4;
%     logic_type(2) = 4;
%     logic_type(3) = 4;
%     logic_type(4) = 4;
%     J(1,2) = 1;
%     J(2,1) = 1;
%     J(3,4) = 1;
%     J(4,3) = 1;

% 6) INDEPENDENT SELF-COPY
%     logic_type(1) = 4;
%     logic_type(2) = 4;
%     logic_type(3) = 4;
%     logic_type(4) = 4;
%     J(1,1) = 1;
%     J(2,2) = 1;
%     J(3,3) = 1;
%     J(4,4) = 1;

    
% ---------------------------------------------------------------

% INDEPENDENT SELF-COPY
%     logic_type(1) = 4;
%     logic_type(2) = 4;
%     logic_type(3) = 4;
%     logic_type(4) = 4;
%     J(1,1) = 1;
%     J(2,2) = 1;
%     J(3,3) = 1;
%     J(4,4) = 1;
    
% IND-NOISY-CPY(2) JOINED BY NOISY AND
%     logic_type(1) = 4;
%     logic_type(2) = 4;
%     logic_type(3) = 1;
%     J(1,1) = 1;
%     J(2,2) = 1;
%     J(3,[1 2]) = 1;

% PARITY AND MINORITY FROM 4 OUTPUT NODES
%     logic_type(1) = 6;
%     logic_type(2) = 6;
%     logic_type(3) = 6;
%     logic_type(4) = 6;
%     logic_type(5) = 9;
%     logic_type(6) = 8;
% 
%     J(5,[1 2 3 4]) = 1;
%     J(6,[1 2 3 4]) = 1;

% MAJORITY
%     logic_type(1) = 7;
%     logic_type(2) = 7;
%     logic_type(3) = 7;
%     logic_type(4) = 7;
%     logic_type(5) = 7;
%     J(1,[1 2 3 4 5]) = 1;
%     J(2,[1 2 3 4 5]) = 1;
%     J(3,[1 2 3 4 5]) = 1;
%     J(4,[1 2 3 4 5]) = 1;
%     J(5,[1 2 3 4 5]) = 1;
    
% MAJORITY - NOT REDUNDANT
%     logic_type(1) = 7;
%     logic_type(2) = 7;
%     logic_type(3) = 7;
%     logic_type(4) = 7;
%     logic_type(5) = 7;
%     J(1,[1 2 3]) = 1;
%     J(2,[1 2 5]) = 1;
%     J(3,[3 4 5]) = 1;
%     J(4,[1 2 4]) = 1;
%     J(5,[2 3 4]) = 1;

% EI VS. MIP TEST
%     logic_type(1) = 4;
%     logic_type(2) = 4;
%     logic_type(3) = 1;
%     logic_type(4) = 4;
%     logic_type(5) = 4;
%     logic_type(6) = 6;
%     logic_type(7) = 6;
% %     logic_type(8) = 6;
% %     logic_type(9) = 6;
% %     logic_type(10) = 4;
% %     logic_type(11) = 4;
% %     logic_type(12) = 4;
% %     logic_type(13) = 4;
% %     logic_type(14) = 4;
%     J(1,[6]) = 1;
%     J(2,[7]) = 1;
%     J(3,[1 2]) = 1;
% %     J(4,[1]) = 1;
%     J(4,[1]) = 1;
%     J(5,[2]) = 1;
% %     J(7,[2]) = 1;
%     J(6,[]) = 1;
%     J(7,[]) = 1;
%     J(10,[4]) = 1;
%     J(11,[5]) = 1;
%     J(12,[3]) = 1;
%     J(13,[6]) = 1;
%     J(14,[7]) = 1;

    
% OPTIMIZED AND GATES (BALDUZZI TONONI 08)
%     logic_type(1) = 1;
%     logic_type(2) = 1;
%     logic_type(3) = 1;
%     logic_type(4) = 1;
%     logic_type(5) = 1;
%     logic_type(6) = 1;
%     logic_type(7) = 1;
%     logic_type(8) = 1;
%     J(1,[2 8]) = 1;
%     J(2,[3 5]) = 1;
%     J(3,[]) = 1;
%     J(4,[6 8]) = 1;
%     J(5,[4 8]) = 1;
%     J(6,[5 7]) = 1;
%     J(7,[1 3]) = 1;
%     J(8,[]) = 1;

% 1 XOR, 2 NOT
%     logic_type(1) = 2;
%     logic_type(2) = 4;
%     logic_type(3) = 3;
%     J(1, [2 3]) = 1;
%     J(2, [3]) = 1;
%     J(3,[1 2]) = 1;

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

% THIS IS NEVER USED SO IT IS NOW COMMENTED OUT SINCE IT'S A HUGE LOOP
% p1 = ones(2^N,2^N); % p(x0,x1) = p(x1|x0)
% 
% for i=1: 2^N
%     for j=1: 2^N
%         x1 = trans2(j-1,N);
%         for k=1: N
%             if x1(k) == 1
%                 p1(i,j) = p1(i,j)*p_x0(i,k);
%             else
%                 p1(i,j) = p1(i,j)*(1-p_x0(i,k));
%             end
%         end
%     end
% end

% p = p1; % practical version
if op_TPM == 0
    p = p_x0; % original version
else
    load TPM;
    N = size(p,2);
end

tpm = p;
save('tpm.mat','tpm')

% if op_disp > 1
%     figure(100)
%     subplot(1,2,1),imagesc(J)
%     colormap('gray')
%     subplot(1,2,2),imagesc(p)
% end
% pause;

% f = zeros(N,1); % averaged firing rates
% for i=1: N
%     f(i) = sum(p_x0(:,i))/2^N;
% end
% avef = sum(f)/N;
% fprintf('avef=%f\n',avef);

if (op_ave == 1)
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
end

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
%     fprintf('x1=%s\n',mat2str(x1));
    
    % partial_prob_comp(partition, partition, state, prob_matrix, binary
    % table, op_fb
    check_prob = partial_prob_comp(1:N,1:N,x1,p,b_table,1); % last argument is op_fb = 1;
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
                if op_context == 0
%                     [BRs FRs] = comp_pers(x1,p,b_table,options);
                    [Big_phi phi prob_cell MIPs M_IRR] = big_phi_comp_fb(M,x1,p,b_table,options);
                    % irreducible points
                    [IRR_REP IRR_phi IRR_MIP M_IRR] = IRR_points(prob_cell,phi,MIPs,M, 0,op_fb);
                    fprintf('\n')
                    fprintf('---------------------------------------------------------------------\n\n')
                    fprintf('Big_phi = %f\n', Big_phi);
                    fprintf('Sum of small_phis = %f\n',sum(phi));
                    fprintf('\nCore Concepts For Complex (Purview, MIP(past & future), Small phi):\n\n');
                    plot_REP(Big_phi, IRR_REP,IRR_phi,IRR_MIP, 1, M, options)
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
        Big_phi_ave = sum(p_x1 .* Big_phi_st); %weighted ave/expected value
    end
    fprintf('Big_phi_ave=%f\n',Big_phi_ave);
end

op_close = 0;
isOpen = matlabpool('size');
if isOpen > 0 && op_close == 1
    matlabpool close;
end

fprintf('\n');

toc;