function output = part_tester(adj_matrix)

nVec = 1:size(adj_matrix,1);
cliques = logical(maximalCliques(adj_matrix));
max_clique_sizes = sum(cliques,1);
max_clique_size_total = max(max_clique_sizes);

clique_sets = cell(max_clique_size_total,1);


for i = 1:size(cliques,2)
    
    for j = 2:max_clique_sizes(i)
        
        concepts = nVec(cliques(:,i));
        clique_sets{j} = union_edit(clique_sets{j}, nchoosek(concepts,j), 'rows');
        
    end
end

output = clique_sets;