function denom_conditional_joint = comp_pers_cpt(num_nodes_indices,denom_nodes_indices,numerator_state,bf_option)

%  compute BRs and FRs for a single perspective but given some fixed
%  current state

if isempty(denom_nodes_indices)
    denom_conditional_joint = [];
    return
% elseif isempty(num_nodes_indices)
% %     num_sys_nodes = denom_nodes_indices(1).num_sys_nodes;
% %     denom_conditional_joint_size = ones(1,2*num_sys_nodes);
% %     denom_conditional_joint_size(1:num_sys_nodes == denom_nodes_indices
%     denom_conditional_joint = [];
%     return
end

global nodes

numerator_nodes = nodes(num_nodes_indices);
num_sys_nodes = nodes(1).num_sys_nodes;
denom_nodes_shift = denom_nodes_indices + num_sys_nodes;
denom_nodes = nodes(denom_nodes_shift);


if strcmp(bf_option,'backward')
    
    
% P(denom_nodes_f | num_nodes_c = numerator_state) = P(denom_nodes_c | num_nodes_p = numerator_state)
elseif strcmp(bf_option,'forward')
    
    % we do the first iteration outside the for main foor loop so we can
    % initialize the joint
    denom_conditional_joint = denom_nodes(1).cpt;
    
    conditioning_indices = cell(1,2*num_sys_nodes);
    for i = 1:2*num_sys_nodes
        conditioning_indices{i} = ':';
    end
    
    % marginalize over nodes not in numerator, these nodes are outside the
    % system for this iteration or they are outside a partition - either
    % way we apply maxent prior/marginalization
    for j = 1:num_sys_nodes
        
        if ~any(j == num_nodes_indices) && any(j == denom_nodes(1).input_nodes)
            denom_conditional_joint = ...
                sum(denom_conditional_joint,j)./size(denom_conditional_joint,j);
        elseif any(j == num_nodes_indices)
            conditioning_indices{j} = numerator_state(j) + 1;
        end
        
    end
    
    for num_i = 2:length(denom_nodes)
        
        next_denom_node_distribution = denom_nodes(num_i).cpt;
        
        % marginalize over nodes not in denom, these nodes are outside the
        % system for this iteration or they are outside a partition - either
        % way we apply maxent prior/marginalization
        for j = 1:num_sys_nodes

            if ~any(j == num_nodes_indices) && any(j == denom_nodes(num_i).input_nodes)
                next_denom_node_distribution = ...
                    sum(next_denom_node_distribution,j)./size(next_denom_node_distribution,j);
            end
        end
        
        % the magic
        denom_conditional_joint = bsxfun(@times,denom_conditional_joint,next_denom_node_distribution);
    end
    
    
    % can we pick these out before hand, or no?
    denom_conditional_joint = denom_conditional_joint(conditioning_indices{:});
    permute_order = [num_sys_nodes+1:2*num_sys_nodes 1:num_sys_nodes];
    denom_conditional_joint = permute(denom_conditional_joint,permute_order);

end







% N = size(p,2);
% 
% % x0 = find(b_table{current,N} == 1);
% x0 = current;
% 
% % x0 = M_pb{current}; % current subset
% x0_out = pick_rest(1:N,x0); % complement
% x0_si = x0_s(x0); % get the current state of the current subset
% 
% if (bf_option == 1) % backwards
% %     xp = M_pb{past}; % past
% %     xp = find(b_table{other,N} == 1);
%     xp = other;
%     xp_out = pick_rest(1:N,xp);
%     rep = partial_prob_cons(xp,[],xp_out,x0,x0_si,p,b_table,1); % all out except for xp
%     
%     if ~isempty(xp)
%         matrix_size = ones(1,N);
%         matrix_size(xp) = 2;
%         rep = reshape(rep,matrix_size);
%     end
%     
%     
% end
% 
% if (bf_option == 2) % forwards
% %     xf = M_pb{future}; % future
% %     xf = find(b_table{other,N} == 1);
%     xf = other;
%     rep = partial_prob_cons(x0,[],x0_out,xf,x0_si,p,b_table,0); % all out except for x0
%     
%     if ~isempty(xf)
%         matrix_size = ones(1,N);
%         matrix_size(xf) = 2;
%         rep = reshape(rep,matrix_size);
%     end
% end
% 
% 
% 
%     
% 
% end