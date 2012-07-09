function big_phi_mip = big_phi_shift(IRR_whole,concepts_whole_p,concepts_whole_f,phi_whole,IRR_parts,concepts_parts_p,concepts_parts_f,phi_parts,op_big_phi_dist)

nWholeConcepts = length(IRR_whole);
nPartConcepts = length(phi_parts);
partitionedCheck = zeros(nPartConcepts,1);

dist_matrix = gen_dist_matrix(size(concepts_whole_p,1));

distance_sum = 0;
phi_sum = 0;


for i = 1:nWholeConcepts
    
    
    IRR_w = IRR_whole{i};
    % For debuggin take out
    fprintf('Concept From Whole: %s\n',mod_mat2str(IRR_w));
%     concept_w = concepts_whole(:,i);
%     phi_w = phi_whole(i);
    
    partner_found = 0;
    
    for j = 1:nPartConcepts
        
        % we haven't already found this guy's partner
        if (partitionedCheck(j) == 0)
            
            IRR_p = IRR_parts{j};
        
            if (length(IRR_p) == length(IRR_w) && all(ismember(IRR_w,IRR_p)))
                
                partner_found = 1;
                partitionedCheck(j) = 1;
                
                if op_big_phi_dist == 0
                    %add in distances b/w past concepts for this purview
                    past_dist = KLD(concepts_whole_p(:,i), concepts_parts_p(:,j));
                    distance_sum = distance_sum + past_dist;
                    %add in distances b/w future concepts for this purview
                    future_dist = KLD(concepts_whole_f(:,i), concepts_parts_f(:,j));
                    distance_sum = distance_sum + future_dist;
                    
                elseif op_big_phi_dist == 1
                    %add in distances b/w past concepts for this purview
                    past_dist = emd_hat_gd_metric_mex(concepts_whole_p(:,i), concepts_parts_p(:,j),dist_matrix);
                    distance_sum = distance_sum + past_dist;
                    %add in distances b/w future concepts for this purview
                    future_dist = emd_hat_gd_metric_mex(concepts_whole_f(:,i), concepts_parts_f(:,j),dist_matrix);
                    distance_sum = distance_sum + future_dist;
                
                end
                
                %for deubbing, take out
                fprintf('\tDistance to past distribution: %f\n',past_dist);
                fprintf('\tDistance to future distribution: %f\n',future_dist);
                fprintf('\tSmall Phi Diff: %f(abs) %f(whole - part)\n',abs(phi_whole(i) - phi_parts(j)),phi_whole(i) - phi_parts(j));
                %add in phi difference
                phi_sum = phi_sum + abs(phi_whole(i) - phi_parts(j));
                
            end
        end
        
    end
    
    % if we didn't find a partner, just add in the phi value
    if ~partner_found
        fprintf('\tConcept does not exist in partitioned system\n');
        fprintf('\tSmall Phi Contribution: %f\n',phi_whole(i));
        phi_sum = phi_sum + phi_whole(i);

    end
end

% for concepts in the partitioned system which do not exist in the whole,
% add their small_phi in
if (any(partitionedCheck == 0))
    fprintf('\tPartitioned concepts not in whole:\n');
end

for i = 1:nPartConcepts
    
    if (partitionedCheck(i) == 0)
        
        fprintf('\t\t%s: small phi contribution, %f',mod_mat2str(IRR_parts{i}),phi_parts(i));
        phi_sum = phi_sum + phi_parts(i);
        
    end
end

big_phi_mip = distance_sum + phi_sum;
                
                

        

