function [Big_phi_M phi_M prob_M subsets MIP_M M_IRR_M] = big_phi_all(network,whole_sys_state)
% [Big_phi_M phi_M prob_M subsets MIP_M] = big_phi_all(x0_s,p,b_table,options)
% this is a test.
%
% this should show up.
% See also: big_phi_comp

% this should not show up

% compute Big-phi in every possible subset
% x0_s : given data about the current neural state
% p: transition matrix in the whole system
% op: the way how small phi is computed (0: difference of entropy, 1: KLD)
% op_disp: (1: display the results, 0: not)

op_console = network.options(10);


% global b_table
% global BRs, global FRs

N = network.num_nodes; % number of elements in the whole system
nodes_vec = network.full_system;

% subset - build a cell array that contains all of the subsets

% subsets builds arrays that use the actual node numbers as opposed to
% logicals - perhaps we should make one of these that is global as well
subsets = cell(network.num_states-1,1); % subtract one since we don't consider the empty system

for i = 1:network.num_states-1 % don't include empty set, this is why for-loop starts at 2
    
    subsets{i} = nodes_vec(logical(network.b_table{i+1,N}));
    
end



% compute big phi in every possible subset
Big_phi_M = zeros(network.num_states-1,1); % Big_phi for each subset except the empty set
phi_M = cell(network.num_states-1,1);
prob_M = cell(network.num_states-1,2); 
MIP_M = cell(network.num_states-1,1); % the partition that gives Big_phi_MIP for each subset
M_IRR_M = cell(network.num_states-1,1);


parfor sub_index = 1:network.num_states-1 % for all non empty subsets of the system\
    
    this_subset = subsets{sub_index}; % get the subset
    if op_console
        fprintf('--------------------------------------------------------------\n\n')
        fprintf('System = %s\n\n',mod_mat2str(this_subset));
    end
    
    [Big_phi phi prob_cell MIP M_IRR] = big_phi_comp_fb(this_subset,whole_sys_state,network); 

    MIP_M{sub_index} = MIP;

    Big_phi_M(sub_index) = Big_phi; % Big_phi for each subset
    phi_M{sub_index} = phi; % Set of small_phis for each purview of each subset
    M_IRR_M{sub_index} = M_IRR; % numerators of purviews with non-zero/positive phi
    
    % concept distributions
    prob_M(sub_index,:) = prob_cell(:); % first layer is subset, second is purview, third is backward/forward
%     prob_M{sub_index,2} = prob_cell{2}; % same as above but for MIP
    
end


% display

if op_console
    fprintf('\n')
    fprintf('--------------------------------------------------------------\n\n')
    fprintf('Big phi values in subset this_subset:\n\n')
    for sub_index = 1: network.num_states-1
        if (Big_phi_M(sub_index) ~= 0 && ~isnan(Big_phi_M(sub_index)))
            fprintf('this_subset=%s: Big_phi=%f\n',mod_mat2str(subsets{sub_index}),Big_phi_M(sub_index));
        end
    end
end
