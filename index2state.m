function state = index2state(index, state_size_vec)

num_nodes = length(state_size_vec);
state = zeros(num_nodes,1);
index = index - 1; % we do this to compensate for the 0 state
for i = 1:num_nodes
%     multidec_array(i) = value - 2 * floor(value/2);
    state(i) = mod(index,state_size_vec(i));
    index = floor(index/state_size_vec(i));
end