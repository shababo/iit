function [Big_phi_MIP MIP] = MIP_search(M,N,Big_phi_M,M_IRR_M,prob_M, phi_M,op_big_phi, op_sum)

%%
% Find the Big-phi MIP in a subset M
% M: a subset where Big_phi_MIP is computed
% N: number of elements in the whole system
% Big_phi_M: Big_phi values in every subset, M
% prob - distributions for the concept for each purview
% phi - phi values for each purview in prob

%%

global grain;

N_M = length(M);

N_Bp = 0;
for i=1: floor(N_M/2)
    N_Bp = N_Bp + nchoosek(N_M,i);
end

Big_phi_cand = zeros(N_Bp,2);
MIP_cand = cell(N_Bp,1);

Big_phi_w = Big_phi_M(trans_M(M,N));

l = 1;
for i=1: floor(N_M/2)
    C = nchoosek(M,i);
    N_C = size(C,1);
    for j=1: N_C
        M1 = C(j,:);
        M2 = pick_rest(M,M1);
        
        M1_i = trans_M(M1,N);
        M2_i = trans_M(M2,N);
        
        if(op_sum == 1 || op_big_phi == 0)
            
            Big_phi_partition = Big_phi_M(M1_i) + Big_phi_M(M2_i);
                      
        elseif(op_big_phi == 1)
            
            phi = [phi_M{M1_i}' phi_M{M2_i}'];

            if (~all(phi == 0))

                concepts = zeros(2^N_M,sum(phi~=0));
                concept_phis = phi(phi ~= 0);

                z = 1;
                for k = 1:length(phi)

                    if (phi(k) ~= 0)
                        
                        if(z <= length(M1_IRR))
                            concepts(:,z) = expand_prob(prob_M{M1_i,1}{k}{1},M,M1);
                        else
                            concepts(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i})}{1},M,M2);
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
            phi = [phi_M{M1_i}' phi_M{M2_i}'];
            
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
                            concepts(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i})}{1},M,M2);
                        end
                        z = z + 1;

                    end

                end

                
                Big_phi_partition = big_phi_info(IRRs,concepts,concept_phis);
                
            else

                Big_phi_partition = 0;
            end
        
        elseif(op_big_phi == 3)
            
            
            
            phi = [phi_M{M1_i}' phi_M{M2_i}'];
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

                concepts = zeros(2^N_M,nIRR);
                concept_phis = phi(phi ~= 0);

                z = 1;
                for k = 1:length(phi)

                    if (phi(k) ~= 0)
                        
                        if(z <= sum(phi_M{M1_i} ~= 0))
                            concepts(:,z) = expand_prob(prob_M{M1_i,1}{k}{1},M,M1);
                        else
                            concepts(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i})}{1},M,M2);
                        end
                        z = z + 1;

                    end

                end
                display = 0;
                Big_phi_partition = big_phi_spacing(concepts,concept_phis,display);

            elseif (sum(phi ~= 0) == 1)
                Big_phi_partition = phi((phi ~= 0));
            else
                Big_phi_partition = 0;
            end

        end
        
        d_Big_phi = Big_phi_w - Big_phi_partition;
        
            
%         
        Norm = min(length(M1),length(M2)) + min(length(M1),length(M2));
        % Norm = 1; % No normalization
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

if (op_big_phi == 0)
    [min_norm_Big_phi i_phi_min] = min(Big_phi_cand(:,1));
else
    [min_norm_Big_phi i_phi_min] = min(Big_phi_cand(:,2));
end
Big_phi_MIP = Big_phi_cand(i_phi_min,1);
MIP = MIP_cand{i_phi_min};
M2 = pick_rest(M,MIP);

fprintf('M = %s, MIP = %s-%s, Big_phi_MIP = %f\n',mat2str(M),mod_mat2str(MIP),mod_mat2str(M2),Big_phi_MIP);