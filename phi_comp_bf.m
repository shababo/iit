function [phi_MIP prob prob_prod_MIP MIP] = phi_comp_bf(options,M,x0,xp,xf,x0_s,p)
% compute small phi of a given purview...?
% 
% options = the options
% M = a system
% x0 = state of the system??
% xp = 
% xf = 
% x0_s = 
% p = TPM as a 2^N x N matrix
% b_table
% BRs
% FRs


% op_fb = options(1);
% op_phi = options(2);
% op_disp = options(3);
% op_single = options(4);
% op_ex = options(5);
op_context = options(6);
op_whole = options(7);
op_min = options(9);
op_normalize = options(14);
op_small_phi = options(16);

global BRs, global FRs, global b_table
global BRs_check, global FRs_check
global func_time, global inline_time
global cpt_time tpm_time
global nodes

N = max(M);
Np = length(xp);
N0 = length(x0);
Nf = length(xf);

%% unpartitioned transition repertoire


% current_convi = convi(x0); past_convi = convi(xp); future_convi = convi(xf);
current = sum(2.^(x0-1))+1; past = sum(2.^(xp-1))+1; future = sum(2.^(xf-1))+1;

% if (current_convi ~= current) || (past_convi ~= past) || (future_convi ~= future)
%     disp('ERROR!')
% end

% current = x0; past = xp; future = xf;

if isempty(BRs{current,past})
%     tic
    BRs{current,past} = comp_pers_cpt(x0,xp,x0_s,'backward');
% cpt_time = cpt_time + toc;
% disp('old')
% tic
%     BRs_check{current,future} = comp_pers_single(x0,xp,x0_s,p,1);
% tpm_time = tpm_time + toc;
%     BRs{current,past} = comp_pers_single(x0,xp,x0_s,p,1);
    
%     if ~all(BRs{current,past} == BRs_check{current,past})
%         disp('BR CHECK:')
%         disp(x0)
%         disp(xp)
%         disp(BRs{current,past})
%         disp(BRs_check{current,past})
%         disp(BRs{current,past}(:) == BRs_check{current,past}(:))
%     end
end
prob_bw = BRs{current,past};

if isempty(FRs{current,future})

% tic
    FRs{current,future} = comp_pers_cpt(x0,xp,x0_s,'forward');
% cpt_time = cpt_time + toc;
% tic
%     FRs_check{current,future} = comp_pers_single(x0,xp,x0_s,p,2);
% tpm_time = tpm_time + toc;
    
%     disp('new result:')
%     disp(size(FRs{current,future}))
% %     disp(FRs{current,future}(:))
%     disp('check result:')
%     disp(size(FRs_check{current,future}))
%     disp(FRs_check{current,future})

%     if ~all(FRs{current,future} == FRs_check{current,future})
%         disp('FR CHECK:')
%         disp(x0)
%         disp(xp)
%         disp(FRs{current,future})
%         disp(FRs_check{current,future})
%         disp(FRs{current,future}(:) == FRs_check{current,future}(:))
%     end
end
prob_fw = FRs{current,future};


prob = cell(2,1);
prob{1} = prob_bw(:);
prob{2} = prob_fw(:);

%% more than one
if Np ~= 0
    [xp_b1 xp_b2 Np_b] = bipartition(xp,Np); % partition of xp
else
    xp_b1{1} = []; xp_b2{1} = []; Np_b = 1;
end
% if Nf ~= 0
%     [xf_b1 xf_b2 Nf_b] = bipartition(xf,Nf); % partition of xf
% else
%     xf_b1{1} = []; xf_b2{1} = []; Nf_b = 1;
% end
[x0_b1 x0_b2 N0_b] = bipartition(x0,N0,1); % partition of x0




phi_cand = zeros(Np_b,N0_b,2,2);
prob_prod_vec = cell(Np_b,N0_b,2);

for i=1: Np_b % past or future
    xp_1 = xp_b1{i};
    xp_2 = xp_b2{i};

    for j=1: N0_b % present
        x0_1 = x0_b1{j};
        x0_2 = x0_b2{j};
        
        Norm = Normalization(xp_1,xp_2,x0_1,x0_2);
        
%         if (all(x0 == [1 2]) && all(xp == [1 2]))
%             disp('x0_1:')
%             disp(x0_1)
%             disp('x0_2:')
%             disp(x0_2)
%             disp('xp_1:')
%             disp(xp_1)
%             disp('xp_2:')
%             disp(xp_2)
%             disp('Norm:')
%             disp(Norm)
%         end
%         current_1_convi = convi(x0_1);
%         current_2_convi = convi(x0_2);
%         other_1_convi = convi(xp_1);
%         other_2_convi = convi(xp_2);
        
        current_1 = sum(2.^(x0_1-1))+1;
        current_2 = sum(2.^(x0_2-1))+1;
        other_1 = sum(2.^(xp_1-1))+1;
        other_2 = sum(2.^(xp_2-1))+1;
        
%         if (current_1_convi ~= current_1) || (current_2_convi ~= current_2) || (other_1_convi ~= other_1) || (other_2 ~= other_2_convi)
%             disp('ERROR!')
%         end

        
%         current_1 = x0_1;
%         current_2 = x0_2;
%         other_1 = xp_1;
%         other_2 = xp_2;

        for bf=1: 2 % past and future

            if Norm ~= 0

                if bf == 1
                    if isempty(BRs{current_1,other_1})
%                             tic
                            BRs{current_1,other_1} = comp_pers_cpt(x0_1,xp_1,x0_s,'backward');
%                             cpt_time = cpt_time + toc;
%                             % disp('old')
%                             tic
%                                 BRs_check{current_1,other_1} = comp_pers_single(x0_1,xp_1,x0_s,p,1);
%                             tpm_time = tpm_time + toc;
 
%                             if ~all(BRs{current_1,other_1} == BRs_check{current_1,other_1})
%                             disp('BR CHECK:')
%                             disp(x0_1)
%                             disp(xp_1)
%                             disp(BRs{current_1,other_1})
%                             disp(BRs_check{current_1,other_1})
%                             disp(BRs{current_1,other_1}(:) == BRs_check{current_1,other_1}(:))
%                             end
                    end
                    prob_p1 = BRs{current_1,other_1};

                    if isempty(BRs{current_2,other_2})
%                             tic
                            BRs{current_2,other_2} = comp_pers_cpt(x0_2,xp_2,x0_s,'backward');
%                             cpt_time = cpt_time + toc;
%                             % disp('old')
%                             tic
%                                 BRs_check{current_2,other_2} = comp_pers_single(x0_2,xp_2,x0_s,p,1);
%                             tpm_time = tpm_time + toc;
%                         if ~all(BRs{current_2,other_2} == BRs_check{current_2,other_2})
%                             disp('BR CHECK:')
%                             disp(x0_2)
%                             disp(xp_2)
%                             disp(BRs{current_2,other_2})
%                             disp(BRs_check{current_2,other_2})
%                             disp(BRs{current_2,other_2}(:) == BRs_check{current_2,other_2}(:))
%                         end
                    end
                    prob_p2 = BRs{current_2,other_2};

                else

                    if isempty(FRs{current_1,other_1})
%                         tic
                        FRs{current_1,other_1} = comp_pers_cpt(x0_1,xp_1,x0_s,'forward');
%                         cpt_time = cpt_time + toc;
%                         tic
%                         FRs_check{current_1,other_1} = comp_pers_single(x0_1,xp_1,x0_s,p,2);
%                         tpm_time = tpm_time + toc;
%     
%                         if ~all(FRs{current_1,other_1} == FRs_check{current_1,other_1})
%                             disp('FR CHECK:')
%                             disp(x0_1)
%                             disp(xp_1)
%                             disp('new result')
%                             disp(FRs{current_1,other_1})
%                             disp('old result')
%                             disp(FRs_check{current_1,other_1})
%                             disp(FRs{current_1,other_1}(:) == FRs_check{current_1,other_1}(:))
%                         end
                    end
                    prob_p1 = FRs{current_1,other_1};

                    if isempty(FRs{current_2,other_2})
%                         tic
                        FRs{current_2,other_2} = comp_pers_cpt(x0_2,xp_2,x0_s,'forward');
%                         cpt_time = cpt_time + toc;
%                         tic
%                         FRs_check{current_2,other_2} = comp_pers_single(x0_2,xp_2,x0_s,p,2);
%                         tpm_time = tpm_time + toc;
%     
%                         if ~all(FRs{current_2,other_2} == FRs_check{current_2,other_2})    
%                             disp('FR CHECK:')
%                             disp(x0_2)
%                             disp(xp_2)
%                             disp(FRs{current_2,other_2})
%                             disp(FRs_check{current_2,other_2})
%                             disp(FRs{current_2,other_2}(:) == FRs_check{current_2,other_2}(:))
%                         end
                    end
                    prob_p2 = FRs{current_2,other_2};

                end
                
                if exist('prob_prod_comp','file') == 3
%                     tic
                    prob_p = prob_prod_comp(prob_p1(:),prob_p2(:),xp,xp_1,0); % ADDED (:)
%                     func_time = func_time + toc;
                
                else    
%                     tic
                    if isempty(prob_p1)
                        prob_p = prob_p2(:);
                    elseif isempty(prob_p2)
                        prob_p = prob_p1(:);
                    else
    %                     prob_p_test = (expand_prob(prob_p1,xp,xp_1) .* expand_prob(prob_p2,xp,xp_2));
    %                     prob_p_test = prob_p_test(:)/sum(prob_p_test);
    %                     prob_p1_reshape = reshape(prob_p1,

    %                     repmat_vec = ones(1,N);
    %                     repmat_vec(xp_2) = 2;
    %                     prob_p1_rep = repmat(prob_p1,repmat_vec);
    %                     repmat_vec = ones(1,N);
    %                     repmat_vec(xp_1) = 2;
    %                     prob_p2_rep = repmat(prob_p2,repmat_vec) ;                 
    %                     prob_p_test = prob_p1_rep .* prob_p2_rep;
    %                     prob_p_test = prob_p_test(:);
                        prob_p_test = bsxfun(@times,prob_p1,prob_p2);
                        prob_p = prob_p_test(:);

    %                     if ~all(prob_p_test == prob_p_test2(:))
    %                         disp('***NOPE***')
    %                     end

    %                     prob_p_test = prob_p1(:) * prob_p2(:)';
    %                     prob_p_test = prob_p_test(:)
                    end
%                     inline_time = inline_time + toc;
                end
                
                
%                 if ~all(abs(prob_p - prob_p_test) <= 1e-5)
%                     disp('PROB CHECK:')
%                 disp(abs(prob_p - prob_p_test) <= 1e-5)
%                     disp('****************')
%                     disp('****************')
%                     disp('present part 1:')
%                     disp(x0_1)
%                     disp('other part 1:')
%                     disp(xp_1)
%                     disp('present part 2:')
%                     disp(x0_2)
%                     disp('other part 2:')
%                     disp(xp_2)
%                     disp(prob_p1)
%                     disp(prob_p2)
%                     disp(prob_p)
%                     disp(prob_p_test)
%                     disp('****************')
% %                 else
% %                     disp('present part 1:')
% %                     disp(x0_1)
% %                     disp('other part 1:')
% %                     disp(xp_1)
% %                     disp('present part 2:')
% %                     disp(x0_2)
% %                     disp('other part 2:')
% %                     disp(xp_2)
% %                     disp(prob_p1)
% %                     disp(prob_p2)
% %                     disp(prob_p)
% %                     disp(prob_p_test)
%                 end

                if (op_small_phi == 0)
                    phi = KLD(prob{bf},prob_p);
                elseif (op_small_phi == 1)
                    phi = emd_hat_gd_metric_mex(prob{bf},prob_p,gen_dist_matrix(length(prob_p)));
                elseif op_small_phi == 2
                    phi = k_distance(prob{bf},prob_p);
                end
                prob_prod_vec{i,j,bf} = prob_p;
            else
                prob_prod_vec{i,j,bf} = [];
                phi = Inf;
            end

            phi_cand(i,j,bf,1) = phi;
            phi_cand(i,j,bf,2) = phi/Norm;
             
%             if (all(x0 == [1 2]) && all(xp == [1 2])) && (Norm ~= 0)
%                 disp('bf:')
%                 disp(bf)
%                 disp('full dist')
%                 disp(prob{bf})
%                 disp('part dist')
%                 disp(prob_p)
%                 disp('phi:')
%                 disp(phi)
%                 disp('phi_norm')
%                 disp(phi/Norm)
%             end
        end
    end
end

MIP = cell(2,2,2);
phi_MIP = zeros(1,2);
prob_prod_MIP = cell(2,1);
for bf = 1: 2
    [phi_MIP(bf) i j] = min2(phi_cand(:,:,bf,1),phi_cand(:,:,bf,2),op_normalize);
    prob_prod_MIP{bf} = prob_prod_vec{i,j,bf};

    MIP{1,1,bf} = xp_b1{i};
    MIP{2,1,bf} = xp_b2{i};
    MIP{1,2,bf} = x0_b1{j};
    MIP{2,2,bf} = x0_b2{j};
end

end



function Norm = Normalization(xp_1,xp_2,x0_1,x0_2,xf_1,xf_2)

if nargin == 4
    Norm = min(length(x0_1),length(xp_2)) + min(length(x0_2),length(xp_1));
else
    Norm = min(length(x0_1),length(xp_2)) + min(length(x0_2),length(xp_1)) ...
        + min(length(x0_1),length(xf_2)) + min(length(x0_2),length(xf_1));
end

end

function [X_min i_min j_min k_min] = min3(X,X2,op_normalize)
X_min = Inf; % minimum of normalized phi (or unnormalized if op_normalize == 0)
X_min2 = Inf; % minimum of phi
i_min = 1;
j_min = 1;
k_min = 1;

if (op_normalize == 1)
    for i=1: size(X,1)
        for j=1: size(X,2)
            for k=1: size(X,3)
                if X(i,j,k) <= X_min && X2(i,j,k) <= X_min2
                    X_min = X(i,j,k);
                    X_min2 = X2(i,j,k);
                    i_min = i;
                    j_min = j;
                    k_min = k;
                end            
            end
        end
    end
else
    for i=1: size(X,1)
        for j=1: size(X,2)
            for k=1: size(X,3)
                if X2(i,j,k) <= X_min
    %                 X_min = X(i,j,k);
                    X_min = X2(i,j,k);
                    i_min = i;
                    j_min = j;
                    k_min = k;
                end            
            end
        end
    end
end
end


function [phi_min_choice i_min j_min] = min2(phi,phi_norm,op_normalize)
phi_norm_min = Inf; % minimum of normalized phi
phi_min = Inf; % minimum of phi
i_min = 1;
j_min = 1;

if (op_normalize == 1 || op_normalize == 2)
    for i=1: size(phi,1)
        for j=1: size(phi,2)
%             if phi_norm(i,j) <= phi_norm_min && phi(i,j) <= phi_min
            if phi_norm(i,j) <= phi_norm_min
                phi_min = phi(i,j);
                phi_norm_min = phi_norm(i,j);
                i_min = i;
                j_min = j;
            end
        end
    end
else
    for i=1: size(phi,1)
        for j=1: size(phi,2)
            if phi(i,j) <= phi_min
                phi_min = phi(i,j);
                phi_norm_min = phi_norm(i,j);
                i_min = i;
                j_min = j;
            end
        end
    end
end

if (op_normalize == 0 || op_normalize == 1)
    phi_min_choice = phi_min;
else
    phi_min_choice = phi_norm_min;
end

end
