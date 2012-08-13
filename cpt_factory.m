function cpt = cpt_factory(this_node, inputs, num_total_nodes, output_noise)

dim_sizes = ones(1,num_total_nodes);
dim_sizes(this_node.num) = this_node.num_states;

dim_sizes([inputs.num]) = [inputs.num_states];

cpt = zeros(dim_sizes);

for i = 1:prod([inputs.num_states])
    
    input_state = dec2multibase(i-1,[inputs.num_states]);
    index_vec = ones(1,num_total_nodes);
    index_vec([inputs.num]) = index_vec([inputs.num]) + input_state';
    
    prob_this_node_on = logic_gates(input_state,this_node.logic_type, output_noise);
    index_cell = num2cell(index_vec);
    cpt(index_cell{:}) = 1 - prob_this_node_on;
    
    index_vec(this_node.num) = 2;
    index_cell = num2cell(index_vec);
    cpt(index_cell{:}) = prob_this_node_on;
end