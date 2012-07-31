function iit_run(tpm, in_J, current_state, in_noise, options)
% IIT_RUN Computes concept, small and big phi, and partition information
% for all subsets of a system (exluding the empty set) over a binary network.
%
%   IIT_RUN(TPM, J, CURRENT_STATE, NOISE, OPTIONS) takes a TPM in
%   state x node form, that is TPM(i,j) is the probability node_i = 1 at 
%   time t+1 given that the system is in state j at time t. J is the
%   connectivity matrix of the network such that J(i,j) = 1 when j has a
%   directed edge to i, and J(i,j) = 0 otherwise. current_state is the
%   state of the system at time t (only used if the options are not set to
%   compute over all states). NOISE is a global noise put on all
%   outgoing messages and must take a value on the interval [0,.5]. OPTIONS
%   is a structure of options for the algoirthm created using the
%   set_options function
%
%   see also set_options

tic

fprintf('\nRunning...\n\n')

% get num_nodes, the number of nodes in the whole system
num_nodes = size(tpm,2);

global grain, global noise, global BRs, global FRs, global J, global b_table
global output_data

grain = 50;
noise = in_noise;
J = in_J;
output_data.tpm = tpm;
output_data.J = J;
output_data.current_state = current_state;
output_data.noise = noise;
output_data.options = options;
output_data.num_nodes = num_nodes;


num_subsets = 2^num_nodes;

% binary table and states list
% from now on, all loops over subsets/bipartitions will use
% b_table for their ordering
b_table = cell(num_subsets,num_nodes);
states = zeros(num_nodes,num_subsets);
for i = 1:num_nodes
    for j = 1:2^i
        b_table{j,i} = trans2(j-1,i); % CONSIDER FLIPPING THIS LR
        if i == num_nodes
            states(:,j) = trans2(j-1,i);
        end
    end
end

% determine if we are computing over all states or just one
op_ave = options(18);
if op_ave == 0
    state_max = 1;
    states(:,1) = current_state;
else
    state_max = num_subsets;
end

% we should deal with different arguments not being included
% if nargin == 4 
%     J = ones(num_nodes);
% elseif nargin == 5
%     J = in_J;
% end


output_data.states = states;

% parallel computing
if matlabpool('size')
    matlabpool close force;
end

op_parallel = options(19);

if op_parallel
    matlabpool;
end

% find main complex (do system partitions)
op_complex = options(15);
% forward-backward option, 07.31.2012 - only use option 3 - do both
op_fb = options(1);
% conservative or progressive - only ever use conservative
op_context = options(6);

% init cell arrays for results
Big_phi_M_st = cell(state_max,1);
Big_phi_MIP_st = cell(state_max,1);
MIP_st = cell(state_max,1);
Complex_st = cell(state_max,1);
prob_M_st = cell(state_max,1);
phi_M_st = cell(state_max,1);
concept_MIP_M_st = cell(state_max,1);
complex_MIP_M_st = cell(state_max,1);
Big_phi_MIP_all_M_st = cell(state_max,1);
complex_MIP_all_M_st = cell(state_maxs,1);
purviews_M_st = cell(state_max,1);


% for each state
for z = 1:state_max
    
    this_state = states(:,z);
    
    % reinit backward rep and forward reps for each state
    BRs = cell(num_subsets); % backward repertoire
    FRs = cell(num_subsets); % forward repertoire

    
    fprintf(['State: ' num2str(this_state') '\n'])
%     fprintf('this_state=%s\n',mat2str(this_state));
    
    % partial_prob_comp(partition, partition, state, prob_matrix, binary
    % table, op_fb
    check_prob = partial_prob_comp(1:num_nodes,1:num_nodes,this_state,tpm,b_table,1); % last argument is op_fb = 1;
    state_check = sum(check_prob);
    if state_check == 0
        fprintf('\tThis state cannot be realized...\n')
        Big_phi_M_st{z} = NaN;
        Big_phi_MIP_st{z} = NaN;
    else
        if op_complex == 0 % only consider whole system
            if op_fb == 2
                options(1) = 0; [Big_phi_f phi_f prob_cell_f] = big_phi_comp(1:num_nodes,this_state,tpm,b_table,options);
                options(1) = 1; [Big_phi_b phi_b prob_cell_b] = big_phi_comp(1:num_nodes,this_state,tpm,b_table,options);
                Big_phi = Big_phi_f + Big_phi_b;
            elseif op_fb == 3 % THIS IS THE ONLY ONE WE DO NOW? BOTH FORWARD AND BACKWARD SIMULTANEOUSLY
                M = 1:num_nodes;
                if op_context == 0
%                     [BRs FRs] = comp_pers(this_state,tpm,b_table,options);
                    [Big_phi phi prob_cell MIPs M_IRR] = big_phi_comp_fb(M,this_state,tpm,b_table,options);
                    % irreducible points
                    [IRR_REP IRR_phi IRR_MIP M_IRR] = IRR_points(prob_cell,phi,MIPs,M, 0,op_fb);
                    fprintf('\n')
                    fprintf('---------------------------------------------------------------------\n\n')
                    fprintf('Big_phi = %f\n', Big_phi);
                    fprintf('Sum of small_phis = %f\n',sum(phi));
                    fprintf('\nCore Concepts For Complex (Purview, MIP(past & future), Small phi):\n\n');
                    plot_REP(Big_phi, IRR_REP,IRR_phi,IRR_MIP, 1, M, options)
                else
                    [Big_phi phi prob_cell MIP prob_cell2] = big_phi_comp_fb(M,this_state,tpm,b_table,options);
                end
            else
                [Big_phi phi prob_cell] = big_phi_comp(1:num_nodes,this_state,tpm,b_table,options);
            end
            Big_phi_st(z) = Big_phi;
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % THE CURRENT SETTINGS TAKE US HERE    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else % find the complex
            fprintf('\tComputing state...\n')
            [MIP Complex Big_phi_M Big_phi_MIP_M prob_M phi_M concept_MIP_M complex_MIP_M M_cell Big_phi_MIP_all_M complex_MIP_M_all purviews_M] ...
                = big_phi_complex(this_state,tpm,options);
            
            if op_fb == 2
%                 % subindex b means backward and f means forward
%                 IRR_phi_b = IRR_phi{1};
%                 IRR_phi_f = IRR_phi{2};
%                 IRR_REP_b = IRR_REP{1};
%                 IRR_REP_f = IRR_REP{2};
%                 M_IRR_b = M_IRR{1};
%                 M_IRR_f = M_IRR{2};
%                 % state dependent big phi and big phi MIP
%                 Big_phi_st(z) = sum(IRR_phi_b)+sum(IRR_phi_f);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % THE CURRENT SETTINGS TAKE US HERE    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                Big_phi_M_st{z} = Big_phi_M;
                Big_phi_MIP_st{z} = Big_phi_MIP_M;
                MIP_st{z} = MIP;
                Complex_st{z} = Complex;
                prob_M_st{z} = prob_M;
                phi_M_st{z} = phi_M;
                concept_MIP_M_st{z} = concept_MIP_M;
                complex_MIP_M_st{z} = complex_MIP_M;
                Big_phi_MIP_all_M_st{z} = Big_phi_MIP_all_M;
                complex_MIP_all_M_st{z} = complex_MIP_M_all;
                purviews_M_st{z} = purviews_M;
                
            end
        end
    end
    
    % pause;
end

% load up output_data struct
output_data.Big_phi_M = Big_phi_M_st;
output_data.Big_phi_MIP = Big_phi_MIP_st;
% KILL THIS ONE BELOW
output_data.MIP = MIP_st;
output_data.Complex = Complex_st;
output_data.concepts_M = prob_M_st;
output_data.small_phi_M = phi_M_st;
output_data.concept_MIP_M = concept_MIP_M_st;
output_data.complex_MIP_M = complex_MIP_M_st;
output_data.M_cell = M_cell;
output_data.Big_phi_MIP_all_M = Big_phi_MIP_all_M_st;
output_data.complex_MIP_all_M = complex_MIP_all_M_st;
output_data.purviews_M = purviews_M_st;


% if op_ave == 1
%     if op_fb == 0
%         Big_phi_ave = sum(Big_phi_st)/2^num_nodes;
%     else
%         Big_phi_ave = sum(p_x1 .* Big_phi_st); %weighted ave/expected value
%     end
% %     fprintf('Big_phi_ave=%f\n',Big_phi_ave);
% end

op_close = 0;
isOpen = matlabpool('size');
if isOpen > 0 && op_close == 1
    matlabpool close;
end

toc

fprintf('Loading GUI... \n');

save('output_data_sample.mat','output_data');

iit_explorer(output_data)