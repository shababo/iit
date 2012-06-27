function [Big_phi_MIP MIP] = MIP_search(M,N,Big_phi_M,prob_M, phi_M,op_volume)

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
        
        if(op_volume == 1)
            phi = [phi_M{M1_i}; phi_M{M2_i}];

            if (~all(phi == 0))

                concepts = zeros(2^N_M,sum(phi ~= 0));
                concept_phis = zeros(1,sum(phi ~= 0));

                for k = 1:2^N_M-1

                    if (phi(k) ~= 0)
                        if(k <= sum(phi_M{M1_i} ~= 0))
                            concepts(:,k) = expand_prob(prob_M{M1_i,1}{k}{1},M,M1);
                        else
                            concepts(:,k) = expand_prob(prob_M{M1_i,1}{k}{1},M,M2);
                        end
                        concept_phis(i) = phi(i);

                    end

                end

            %     disp(concepts);
            %     disp(concept_phis);
            %     disp(phi);

                Big_phi_partition = big_phi_volume(concepts,concept_phis,grain);

            else
                Big_phi_partition = 0;
            end
        else
            Big_phi_partition = Big_phi_M(M1_i) + Big_phi_M(M2_i);
        end
        

%         
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

if (op_volume == 0)
    [min_norm_Big_phi i_phi_min] = min(Big_phi_cand(:,1));
else
    [min_norm_Big_phi i_phi_min] = min(Big_phi_cand(:,2));
end
Big_phi_MIP = Big_phi_cand(i_phi_min,1);
MIP = MIP_cand{i_phi_min};
M2 = pick_rest(M,MIP);

fprintf('M = %s, MIP = %s-%s, Big_phi_MIP = %f\n',mat2str(M),mod_mat2str(MIP),mod_mat2str(M2),Big_phi_MIP);