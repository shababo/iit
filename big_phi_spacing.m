function big_phi = big_phi_spacing(concepts,c_phi)

nConcepts = size(concepts,2);
dim = size(concepts,1);

emd_dist_matrix = gen_dist_matrix(dim);

overlap = 0;

% compute pairwise distances and if less than .5, subtract appropriate
% amount
for i = 1:nConcepts
    for j = i+1:nConcepts

        dist = emd_hat_gd_metric_mex(concepts(:,i),concepts(:,j),emd_dist_matrix);
       
        if (dist < .5)
            
            overlap = overlap + (.5 - dist) * min(c_phi(i),c_phi(j));
            
        end
    end
end

scaling_factor = nConcepts / nchoosek(nConcepts,2);

big_phi = nConcepts - overlap * scaling_factor;

end