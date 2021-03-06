function [Big_phi_MIP MIP Big_phi_cand MIP_cand] = MIP_search_reentry(M,N,Big_phi_M,M_IRR_M,prob_M, phi_M,options, concept_MIP_M, network)

%%
% Find the Big-phi MIP in a subset M
% M: a subset where Big_phi_MIP is computed
% N: number of elements in the whole system
% Big_phi_M: Big_phi values in every subset, M
% prob - distributions for the concept for each purview
% phi - phi values for each purview in prob

%%

% save('prob_M.mat','prob_M')

% global grain;
% %debug remove me
% fprintf('----------------------------------------------\n');
% fprintf('M = %s\n',mat2str(M));

op_big_phi = options(11);
op_sum = options(12);
op_normalize = options(13);
op_big_phi_dist = options(17);
op_console = options(10);

N_M = length(M);

N_Bp = 0;
for i=1: floor(N_M/2)
    N_Bp = N_Bp + nchoosek(N_M,i);
end

Big_phi_cand = zeros(N_Bp,2);
MIP_cand = cell(N_Bp,1);

whole_i = trans_M(M,N);
Big_phi_w = Big_phi_M(whole_i);

if (any(op_big_phi == [0 6 7 4 5])) %Larissa: For 0 I need to know the concepts it has, but actually not the distributions...
    
    phi_whole = phi_M{whole_i}(:,1)';
    concept_numind = find(phi_whole ~= 0);
    phi_w_concepts = phi_whole(phi_whole ~= 0);
    IRR_whole = M_IRR_M{whole_i};
    if any(op_big_phi == [4 5])
        concepts_whole_p = zeros(2^N_M,length(phi_w_concepts));
        concepts_whole_f = zeros(2^N_M,length(phi_w_concepts));

        z = 1;
        for i = 1:length(phi_whole)  %Larissa: Why not length(phi_w_concepts)?
            if (phi_whole(i) ~= 0)

                if ~isempty(prob_M{whole_i,1}{i}{1})
                    concepts_whole_p(:,z) = prob_M{whole_i,1}{i}{1};
                end
                if ~isempty(prob_M{whole_i,1}{i}{2})
                    concepts_whole_f(:,z) = prob_M{whole_i,1}{i}{2};
                end
                z = z + 1;
            end
        end  

        whole_concepts_cell_p = cell(length(phi_w_concepts),2);
        whole_concepts_cell_f = cell(length(phi_w_concepts),2);

        for i = 1:length(phi_w_concepts)

           whole_concepts_cell_p{i,1} = concepts_whole_p(:,i);  %Larissa: do I need these?
           whole_concepts_cell_f{i,1} = concepts_whole_f(:,i); 
           whole_concepts_cell_p{i,2} = phi_w_concepts(i);
           whole_concepts_cell_f{i,2} = phi_w_concepts(i);

        end
    end
%     
%     phi_whole = phi_w_concepts;
    
end
    

% are we doing both sides of the partition!?!? <-- YES WHEN IT IS HALF AND
% HALF
l = 1;
for i=1: floor(N_M/2)
    C = nchoosek(M,i);
    N_C = size(C,1);
    for j=1: N_C
        M1 = C(j,:);
        M2 = pick_rest(M,M1);
        
        M1_i = trans_M(M1,N);
        M2_i = trans_M(M2,N);
        
        %debug remove me
%         fprintf('Partition: %s x %s\n',mod_mat2str(M1),mod_mat2str(M2)); 
        
        if(op_sum == 1 || op_big_phi == 0 || op_big_phi == 6 || op_big_phi == 7)
            
            %Big_phi_partition = Big_phi_M(M1_i) + Big_phi_M(M2_i);
            PhiCutSum = [0; 0];  %Larissa: cutting first M1 <- M2 (causes on M1, effects from M2) and then M1 -> M2 (causes on M2, effects from M1)
            if (any(op_big_phi == [6 7]))
                BRcut_dist = cell(length(phi_w_concepts), 2, 2); %dim1: per concept, dim2: past/future, dim2: whole/cut
                BRcut_phi = zeros(length(phi_w_concepts),1);
                FRcut_dist = cell(length(phi_w_concepts), 2, 2);
                FRcut_phi = zeros(length(phi_w_concepts),1);
            end    
            for k = 1:length(phi_w_concepts)
                IRR_w = IRR_whole{k};
                if all(ismember(IRR_w,M1)) 
                    % for M1 <- M2 cut take BR of M1 and FR from M
                    
                    %Larissa: As long as the concepts are still not
                    %ordered according to trans_M...
                    % Now definitely NOT NICE!
                    indm = 0;
                    m = 1;
                    mfound = 0;
                    for subset_size = 1:length(M1)
                        M1_purviews = nchoosek(M1,subset_size);
                        N_M1_pur = size(M1_purviews,1);
                        for q = 1:N_M1_pur % for all combos of size i
                            if numel(M1_purviews(q,:)) == numel(IRR_w) && all(M1_purviews(q,:) == IRR_w)
                               indm = m; 
                               break;
                            end   
                            m = m + 1;
                        end
                        if mfound, break, end;
                    end
                                        
                    phi_BRcut = min(phi_M{M1_i}(indm,2), phi_M{whole_i}(concept_numind(k),3));
                    phi_FRcut = min(phi_M{whole_i}(concept_numind(k),2), phi_M{M1_i}(indm,3));
                    if (any(op_big_phi == [6 7])) %L1 or Earthmover
                        %Larissa: distributions that are identical anyways are empty
%                         denom_p = sort([concept_MIP_M{M1_i}{indm}{:,1,1}]);
%                         denom_f = sort([concept_MIP_M{M1_i}{indm}{:,1,2}]);
                        
                        cutpdist = expand_prob(prob_M{M1_i,1}{indm}{1},M,M1);
                        BRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} cutpdist};
                                               
                        cutfdist = expand_prob(prob_M{M1_i,1}{indm}{2},M,M1);   %Larissa: Check Expand Probabilities for future. I'm still not convinced this is correct
                        % The above is not correct, the result here should
                        % not be a flat max ent but a forward maxent...
                        FRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} cutfdist};
                        BRcut_phi(k) = phi_BRcut;
                        FRcut_phi(k) = phi_FRcut;
                        if op_big_phi == 7
                            %Larissa: Maybe more efficient to leave it [] but hopefully doesn't matter much
                            BRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} prob_M{whole_i,1}{concept_numind(k)}{2}}; 
                            FRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} prob_M{whole_i,1}{concept_numind(k)}{1}}; 
                        end    
                    end    
                elseif all(ismember(IRR_w,M2))
                    
                    indm = 0;
                    m = 1;
                    mfound = 0;
                    for subset_size = 1:length(M2)
                        M2_purviews = nchoosek(M2,subset_size);
                        N_M2_pur = size(M2_purviews,1);
                        for q = 1:N_M2_pur % for all combos of size i
                            if numel(M2_purviews(q,:)) == numel(IRR_w) && all(M2_purviews(q,:) == IRR_w)
                               indm = m; 
                               mfound = 1;
                               break;
                            end   
                            m = m + 1;
                        end
                        if mfound, break, end;
                    end
                    
                    phi_BRcut = min(phi_M{whole_i}(concept_numind(k),2), phi_M{M2_i}(indm,3));
                    phi_FRcut = min(phi_M{M2_i}(indm,2), phi_M{whole_i}(concept_numind(k),3));
                    if (any(op_big_phi == [6 7])) %L1 or Earthmover
                        %Larissa: distributions that are identical anyways
                        %are empty for option 6
%                         denom_p = sort([concept_MIP_M{M2_i}{indm}{:,1,1}]);
%                         denom_f = sort([concept_MIP_M{M2_i}{indm}{:,1,2}]);
                        cutfdist = expand_prob(prob_M{M2_i,1}{indm}{2},M,M2);   %wrong expand_prob
                        BRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} cutfdist}; %future might have changed

                        cutpdist = expand_prob(prob_M{M2_i,1}{indm}{1},M,M2); 
                        FRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} cutpdist}; %back might have changed, future is the same
                        BRcut_phi(k) = phi_BRcut;
                        FRcut_phi(k) = phi_FRcut;
                        if op_big_phi == 7
                            %Larissa: Maybe more efficient to leave it [] but hopefully doesn't matter much
                            BRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} prob_M{whole_i,1}{concept_numind(k)}{1}}; %back is the same
                            FRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} prob_M{whole_i,1}{concept_numind(k)}{2}}; %back might have changed, future is the same
                       end 
                    end   
                else % if numerator has elements from both sides
                    denom_p = sort([concept_MIP_M{whole_i}{concept_numind(k)}{:,1,1}]);    %Larissa: The sort may be important 
                    denom_f = sort([concept_MIP_M{whole_i}{concept_numind(k)}{:,1,2}]);
                    % for BRcut (M1 <- M2 is cut) M1M2/[M1]p[M2]f is still intact
                    BRcut_pdist = []; BRcut_fdist = []; FRcut_pdist = []; FRcut_fdist = []; %Only needed for op_big_phi = 6 or 7
                    if all(ismember(denom_p,M1))
                        phi_BRcut_BR = phi_M{whole_i}(concept_numind(k),2); %stays the same
                        if all(ismember(denom_f,M2))
                            phi_BRcut_FR = phi_M{whole_i}(concept_numind(k),3);
                        else
                            % check the new forward phi_mip M1M2/[M1M2]f for all
                            % possible denominators
                            [phi_BRcut_FR, BRcut_fdist, denom_fnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'forward','BRcut');
                        end    
                    else
                        if all(ismember(denom_f,M2))
                            phi_BRcut_FR = phi_M{whole_i}(concept_numind(k),3);
                            % check the new backward phi_mip M1M2/[M1M2]p for all
                            % possible denominators
                            [phi_BRcut_BR, BRcut_pdist, denom_pnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'backward','BRcut');
                        else
                            % check the new back and forward phi_mip M1M2/[M1M2]p[M1M2]f for all
                            % possible denominators
                            [phi_BRcut_BR, BRcut_pdist, denom_pnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'backward','BRcut');
                            [phi_BRcut_FR, BRcut_fdist, denom_fnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'forward','BRcut');
                        end    
                    end
                                        
                    % for FRcut (M1 -> M2 is cut) M1M2/[M2]p[M1]f is still intact
                    if all(ismember(denom_p,M2))
                        phi_FRcut_BR = phi_M{whole_i}(concept_numind(k),2); %stays the same
                        if all(ismember(denom_f,M1))
                            phi_FRcut_FR = phi_M{whole_i}(concept_numind(k),3);
                        else
                            % check the new forward phi_mip M1M2/[M1M2]f for all
                            % possible denominators
                            [phi_FRcut_FR, FRcut_fdist, denom_fnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'forward','FRcut');
                        end    
                    else
                        if all(ismember(denom_f,M1))
                            phi_FRcut_FR = phi_M{whole_i}(concept_numind(k),3);
                            % check the new backward phi_mip M1M2/[M1M2]p for all
                            % possible denominators
                            [phi_FRcut_BR, FRcut_pdist, denom_pnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'backward','FRcut');
                        else
                            % check the new back and forward phi_mip M1M2/[M1M2]p[M1M2]f for all
                            % possible denominators
                            [phi_FRcut_BR, FRcut_pdist, denom_pnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'backward','FRcut');
                            [phi_FRcut_FR, FRcut_fdist, denom_fnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'forward','FRcut');
                        end    
                    end
                    
                    %Larissa: if we would want to take op_big_phi into
                    %account still, this has to be changed
                    phi_BRcut = min(phi_BRcut_BR, phi_BRcut_FR);
                    phi_FRcut = min(phi_FRcut_BR, phi_FRcut_FR);
%                     if phi_BRcut ~= 0
%                         [M1; M2; IRR_w; denom_p; denom_pnew]
%                     end
                    if (any(op_big_phi == [6 7])) %L1 or Earthmover
                        if ~isempty(BRcut_pdist)
                            BRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} BRcut_pdist};
                        elseif op_big_phi == 7
                            BRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} prob_M{whole_i,1}{concept_numind(k)}{1}};
                        end
                        if ~isempty(BRcut_fdist)
                            BRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} BRcut_fdist};
                        elseif op_big_phi == 7    
                            BRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} prob_M{whole_i,1}{concept_numind(k)}{2}};
                        end
                        if ~isempty(FRcut_pdist)
                            FRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} FRcut_pdist};
                        elseif op_big_phi == 7    
                            FRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} prob_M{whole_i,1}{concept_numind(k)}{1}};
                        end
                        if ~isempty(FRcut_fdist)
                            FRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} FRcut_fdist};
                        elseif op_big_phi == 7    
                                FRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} prob_M{whole_i,1}{concept_numind(k)}{2}};
                        end
                        
                        BRcut_phi(k) = phi_BRcut;
                        FRcut_phi(k) = phi_FRcut;
                    end   
                end % if 
                PhiCutSum = PhiCutSum + [phi_BRcut; phi_FRcut];
             end %for k 
             Big_phi_partition = max(PhiCutSum);       % max here means the minimum of the difference between whole and partitioned system   
        elseif(op_big_phi == 1)
            
            phi = [phi_M{M1_i}(:,1)' phi_M{M2_i}(:,1)'];

            if (~all(phi == 0))

                concepts = zeros(2^N_M,sum(phi~=0));
                concept_phis = phi(phi ~= 0);

                z = 1;
                for k = 1:length(phi)

                    if (phi(k) ~= 0)
                        
                        if(z <= length(M1_IRR))
                            concepts(:,z) = expand_prob(prob_M{M1_i,1}{k}{1},M,M1); 
                        else
                            concepts(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{1},M,M2); 
                        end
                        z = z + 1;

                    end

                end

%         disp('!!!!!!!!!!!!!');
%         disp('!!!!!!!!!!!!!');
%         disp('!!!!!!!!!!!!!');
%         disp(concepts);
%         disp(concept_phis);
%         disp('!!!!!!!!!!!!!');
%         disp('!!!!!!!!!!!!!');
%         disp('!!!!!!!!!!!!!');

                Big_phi_partition = big_phi_volume(concepts,concept_phis,grain);

            else
                Big_phi_partition = 0;
            end
            
        elseif (op_big_phi == 2)
            
            M1_IRR = M_IRR_M{M1_i};
            M2_IRR = M_IRR_M{M2_i};

            nIRR = length(M1_IRR) + length(M2_IRR);
            IRRs = cell(nIRR,1);
            phi = [phi_M{M1_i}(:,1)' phi_M{M2_i}(:,1)'];
            
            for x = 1:nIRR
                if (x <= length(M1_IRR))
                   IRRs{x} = M1_IRR{x};
                else
                   IRRs{x} = M2_IRR{x - length(M1_IRR)};
                end
            end
            

            if (~all(phi == 0))

                concepts = zeros(2^N_M,nIRR);
                concept_phis = phi(phi ~= 0);

                z = 1;
                for k = 1:length(phi)

                    if (phi(k) ~= 0)
                        
                        if(z <= length(M1_IRR))
                            concepts(:,z) = expand_prob(prob_M{M1_i,1}{k}{1},M,M1);
                        else
                            concepts(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{1},M,M2);
                        end
                        z = z + 1;

                    end

                end

                
                Big_phi_partition = big_phi_info(IRRs,concepts,concept_phis);
                
            else

                Big_phi_partition = 0;
            end
        
        elseif(op_big_phi == 3)
            
            
            
            phi = [phi_M{M1_i}(:,1)' phi_M{M2_i}(:,1)'];
            nIRR = sum(phi ~= 0);
            
%             disp('***************')
%             disp(M1_IRR);
%             disp(length(M1_IRR));
%             disp(M2_IRR);
%             disp(length(M2_IRR));
%             disp(nIRR)
%             disp(phi)
%             disp(phi ~= 0);
%             disp(sum(phi ~= 0));
            
            if (nIRR > 1)

                concepts_past = zeros(2^N_M,nIRR);
                concepts_future = zeros(2^N_M,nIRR);
                concept_phis = phi(phi ~= 0);

                z = 1;
                for k = 1:length(phi)

                    if (phi(k) ~= 0)
                        
                        if(z <= sum(phi_M{M1_i}(:,1) ~= 0))
                            concepts_past(:,z) = expand_prob(prob_M{M1_i,1}{k}{1},M,M1);
                            concepts_future(:,z) = expand_prob(prob_M{M1_i,1}{k}{2},M,M1);
                        else
                            concepts_past(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{1},M,M2);
                            concepts_future(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{2},M,M2);
                        end
                        z = z + 1;

                    end

                end
                display = 0;
                Big_phi_partition = big_phi_spacing(concepts_past,concept_phis,display) + big_phi_spacing(concepts_future,concept_phis,display);

            elseif (sum(phi ~= 0) == 1)
                Big_phi_partition = phi((phi ~= 0));
            else
                Big_phi_partition = 0;
            end
            
        elseif (op_big_phi == 4)
            
            M1_IRR = M_IRR_M{M1_i};
            M2_IRR = M_IRR_M{M2_i};

            nIRR = length(M1_IRR) + length(M2_IRR);
            IRR_parts = cell(nIRR,1);
            phi_parts_all = [phi_M{M1_i}(:,1)' phi_M{M2_i}(:,1)'];
            
            for x = 1:nIRR
                if (x <= length(M1_IRR))
                   IRR_parts{x} = M1_IRR{x};
                else
                   IRR_parts{x} = M2_IRR{x - length(M1_IRR)};
                end
            end
            

                
            concepts_past = zeros(2^N_M,nIRR);
            concepts_future = zeros(2^N_M,nIRR);
            phi_parts = phi_parts_all(phi_parts_all ~= 0);

            z = 1;
            
            for k = 1:length(phi_parts_all)

                if (phi_parts_all(k) ~= 0)

                    if(z <= sum(phi_M{M1_i}(:,1) ~= 0))
                        if ~isempty(prob_M{M1_i,1}{k}{1})
                            concepts_past(:,z) = expand_prob(prob_M{M1_i,1}{k}{1},M,M1);
                        end
                        if ~isempty(prob_M{M1_i,1}{k}{2})
                            concepts_future(:,z) = expand_prob(prob_M{M1_i,1}{k}{2},M,M1);
                        end
                    else
                        if ~isempty(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{1})
                            concepts_past(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{1},M,M2);
                        end
                        if ~isempty(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{2})
                            concepts_future(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{2},M,M2);
                        end
                    end
                    
                    z = z + 1;

                end

            end

%                 concepts_parts = zeros(2^N_M,nIRR);
%                 phi_parts = phi(phi ~= 0);
% 
%                 z = 1;
%                 for k = 1:length(phi)
% 
%                     if (phi(k) ~= 0)
%                         
%                         if(z <= length(M1_IRR))
%                             concepts_parts(:,z) = expand_prob(prob_M{M1_i,1}{k}{1},M,M1);
%                         else
%                             concepts_parts(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i})}{1},M,M2);
%                         end
%                         z = z + 1;
% 
%                     end
% 
%                 end

%             fprintf('-------------------------------------------------------\n');
%             fprintf('M = %s with partition %s - %s\n',mod_mat2str(M),mod_mat2str(M1),mod_mat2str(M2));
            d_Big_phi = big_phi_shift(M1_IRR, M2_IRR, N, M, IRR_whole,concepts_whole_p,concepts_whole_f,phi_w_concepts, M1, M2,...
                              IRR_parts,concepts_past,concepts_future, prob_M, phi_M{M1_i}(:,1)', phi_M{M2_i}(:,1)', phi_parts,op_big_phi_dist);
                          
        elseif(op_big_phi == 5)
            
            
            
            phi = [phi_M{M1_i}(:,1)' phi_M{M2_i}(:,1)'];
            nIRR = sum(phi ~= 0);
            
%             disp('***************')
%             disp(M1_IRR);
%             disp(length(M1_IRR));
%             disp(M2_IRR);
%             disp(length(M2_IRR));
%             disp(nIRR)
%             disp(phi)
%             disp(phi ~= 0);
%             disp(sum(phi ~= 0));
            


            concepts_past = zeros(2^N_M,nIRR);
            concepts_future = zeros(2^N_M,nIRR);
            concept_phis = phi(phi ~= 0);

            z = 1;
            for k = 1:length(phi)

                if (phi(k) ~= 0)

                    if(z <= sum(phi_M{M1_i}(:,1) ~= 0))
                        concepts_past(:,z) = expand_prob(prob_M{M1_i,1}{k}{1},M,M1);
                        concepts_future(:,z) = expand_prob(prob_M{M1_i,1}{k}{2},M,M1);
                    else
                        concepts_past(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{1},M,M2);
                        concepts_future(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{2},M,M2);
                    end
                    z = z + 1;

                end

            end
            
            part_concepts_cell_p = cell(nIRR,2);
            part_concepts_cell_f = cell(nIRR,2);

            for k = 1:nIRR

               part_concepts_cell_p{k,1} = concepts_past(:,k);
               part_concepts_cell_f{k,1} = concepts_future(:,k); 
               part_concepts_cell_p{k,2} = concept_phis(k);
               part_concepts_cell_f{k,2} = concept_phis(k);

            end
            
            d_Big_phi = C_distance(part_concepts_cell_p,whole_concepts_cell_p)...
                + C_distance(part_concepts_cell_f,whole_concepts_cell_f);
        
        end
        
        if (any(op_big_phi == [0 1 2 3]))
            
            % 07-06-12 CHANGED TO ABS VALUE
            d_Big_phi = abs(Big_phi_w - Big_phi_partition);
            
        end
        
        if op_big_phi == 6
            back_maxent = expand_prob([],M,[]);
            forward_maxent = expand_prob([],M,[]); %Larissa: WRONG!!

            BRcut_Phi = 0;
            FRcut_Phi = 0;
            for k = 1:size(BRcut_phi)
                if ~isempty(BRcut_dist{k,1,1}) %backward repertoires (if empty they stayed the same!)
                    BRcut_Phi = BRcut_Phi + L1norm(BRcut_dist{k,1,1},BRcut_dist{k,1,2})*BRcut_phi(k) + L1norm(BRcut_dist{k,1,1},back_maxent)*(phi_w_concepts(k)-BRcut_phi(k));
                end
                if ~isempty(BRcut_dist{k,2,1}) %forward repertoires
                    BRcut_Phi = BRcut_Phi + L1norm(BRcut_dist{k,2,1},BRcut_dist{k,2,2})*BRcut_phi(k) + L1norm(BRcut_dist{k,2,1},forward_maxent)*(phi_w_concepts(k)-BRcut_phi(k));
                end
                
                if ~isempty(FRcut_dist{k,1,1}) %backward repertoires (if empty they stayed the same!)
                    FRcut_Phi = FRcut_Phi + L1norm(FRcut_dist{k,1,1},FRcut_dist{k,1,2})*FRcut_phi(k) + L1norm(FRcut_dist{k,1,1},back_maxent)*(phi_w_concepts(k)-FRcut_phi(k));
                end
                if ~isempty(FRcut_dist{k,2,1}) %forward repertoires
                    FRcut_Phi = FRcut_Phi + L1norm(FRcut_dist{k,2,1},FRcut_dist{k,2,2})*FRcut_phi(k) + L1norm(FRcut_dist{k,2,1},forward_maxent)*(phi_w_concepts(k)-FRcut_phi(k));
                end
            end 
            d_Big_phi = min(BRcut_Phi, FRcut_Phi);
        end    
%         
        if op_big_phi == 7  %earth movers for concepts
            back_maxent = expand_prob([],M,[]);
            forward_maxent = expand_prob([],M,[]); %Larissa: WRONG!!
%             indBR = find(BRcut_phi);
%             indFR = find(FRcut_phi);
%             BRcut_phi = BRcut_phi(indBR);
%             FRcut_phi = FRcut_phi(indFR);
%             BRcut_concepts = [BRcut_phi sum(phi_w_concepts)-sum(BRcut_phi)];
%             
%             BRDistMat_past = EMDDistanceMatrix(BRcut_dist(:,1,1), BRcut_dist(indBR,1,2), back_maxent); %past whole and cut distributions
%             BRDistMat_fut = EMDDistanceMatrix(BRcut_dist(:,2,1), BRcut_dist(indFR,2,2), forward_maxent); %future whole and cut distributions
            TempMp = genEMDDistanceMatrix(BRcut_dist(:,1,1), BRcut_dist(:,1,2), back_maxent); %past whole and cut distributions
            BRDistMat_past = [zeros(size(TempMp)), TempMp; TempMp' zeros(size(TempMp))];
            TempMf = genEMDDistanceMatrix(BRcut_dist(:,2,1), BRcut_dist(:,2,2), forward_maxent); %future whole and cut distributions
            BRDistMat_fut = [zeros(size(TempMf)), TempMf; TempMf' zeros(size(TempMf))];

            BRphiDiff = sum(phi_w_concepts)-sum(BRcut_phi);
            tempVphi = [phi_w_concepts'; 0];
            tempVphicut = [BRcut_phi; BRphiDiff];
            BRcut_Phi(1) = emd_hat_gd_metric_mex([tempVphi; zeros(size(tempVphi))],[zeros(size(tempVphicut));tempVphi],BRDistMat_past);
            BRcut_Phi(2) = emd_hat_gd_metric_mex([tempVphi; zeros(size(tempVphi))],[zeros(size(tempVphicut));tempVphi],BRDistMat_fut);
            %BRcut_Phi
            BRcut_Phi = sum(BRcut_Phi); 
            
            
            TempMp = genEMDDistanceMatrix(FRcut_dist(:,1,1), FRcut_dist(:,1,2), back_maxent); %past whole and cut distributions
            FRDistMat_past = [zeros(size(TempMp)), TempMp; TempMp' zeros(size(TempMp))];
            TempMf = genEMDDistanceMatrix(FRcut_dist(:,2,1), FRcut_dist(:,2,2), forward_maxent); %future whole and cut distributions
            FRDistMat_fut = [zeros(size(TempMf)), TempMf; TempMf' zeros(size(TempMf))];
            
            FRphiDiff = sum(phi_w_concepts)-sum(FRcut_phi);

            tempVphicut = [FRcut_phi; FRphiDiff];
            FRcut_Phi(1) = emd_hat_gd_metric_mex([tempVphi; zeros(size(tempVphi))],[zeros(size(tempVphicut));tempVphi],FRDistMat_past);
            FRcut_Phi(2) = emd_hat_gd_metric_mex([tempVphi; zeros(size(tempVphi))],[zeros(size(tempVphicut));tempVphi],FRDistMat_fut);
            %FRcut_Phi
            FRcut_Phi = sum(FRcut_Phi); 
          
            d_Big_phi = min(BRcut_Phi, FRcut_Phi);
        end 
        
        % Norm = min(length(M1),length(M2)) + min(length(M1),length(M2));
        % Norm = 1; % No normalization
        Norm = 2^min(length(M1),length(M2))-1;
        Big_phi_cand(l,1) = d_Big_phi;
        Big_phi_cand(l,2) = d_Big_phi/Norm;
        MIP_cand{l} = M1;
        
%         if N_M == N
%             fprintf('M1=%s-%s: ',mod_mat2str(M1),mod_mat2str(M2));
%             fprintf('%f-(%f+%f)=%f %f\n',Big_phi_w,Big_phi1,Big_phi2,d_Big_phi,d_Big_phi/Norm);
%         end
        
        l = l + 1;
    end
    
end

% if length(M) == N
%     
%     x = cat(1,concepts_whole_p',concepts_past');
%     nWholeConcepts = size(concepts_whole_p,2);
%     save('sample_partition_4n_sys.mat','x','nWholeConcepts')
% %     figure(1)
% %     conceptscatter(all_concepts,size(concepts_whole_p,2))
%     
% end
    

if (op_normalize == 1 || op_normalize == 2) % Option to normalize or not for new methods of computing big phi
    [min_norm_Big_phi i_phi_min] = min(Big_phi_cand(:,2));
else
    [min_norm_Big_phi i_phi_min] = min(Big_phi_cand(:,1));
end

if (op_normalize == 0 || op_normalize == 1)
    Big_phi_MIP = Big_phi_cand(i_phi_min,1);
else
    Big_phi_MIP = Big_phi_cand(i_phi_min,2);
end

MIP = MIP_cand{i_phi_min};
M2 = pick_rest(M,MIP);

if (op_console && Big_phi_MIP ~= 0)
    fprintf('M = %s\nMIP = %s-%s, Big_phi_MIP = %f\n',mat2str(M),mod_mat2str(MIP),mod_mat2str(M2),Big_phi_MIP);
end