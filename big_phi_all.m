function [Big_phi_M phi_M prob_M M_cell MIP_M] = big_phi_all(x0_s,p,b_table,options)
% [Big_phi_M phi_M prob_M M_cell MIP_M] = big_phi_all(x0_s,p,b_table,options)
% this is a test.
%
% this should show up.
% See also: big_phi_comp

% this should not show up

%% compute Big-phi in every possible subset
% x0_s : given data about the current neural state
% p: transition matrix in the whole system
% op: the way how small phi is computed (0: difference of entropy, 1: KLD)
% op_disp: (1: display the results, 0: not)

op_fb = options(1);
op_context = options(6);

N = log2(size(p,1)); % number of elements in the whole system

%% subset - build a cell array that contains all of the subsets
M_cell = cell(2^N-1,1); % subtract one since we don't consider the empty system

for i=1: 2^N-1
    x = trans2(i,N); % we have b_table for this, no?
    C = [];
    for j=1: N % this can all be done with a any() or something like that
        if x(j) == 1
            C = [C j];
        end
    end
    M_cell{i} = C;
end

%% compute big phi in every possible subset
Big_phi_M = zeros(2^N-1,1); % Big_phi for each subset except the empty set
phi_M = cell(2^N-1,1);
prob_M = cell(2^N-1,2); 
MIP_M = cell(2^N-1,1); % the partition that gives Big_phi_MIP for each subset

%% For only conservative, compute BRs and FRs in every possible perspectives 
if op_context == 0
   [BRs FRs] = comp_pers(x0_s,p,b_table,options);
end

%%

for M_i = 1: 2^N-1 % for all proper subsets of the system
    M = M_cell{M_i}; % get the subset
    fprintf('M=%s\n',mod_mat2str(M));
    if op_fb == 3 % YES - fb simultaneously
        if op_context == 0 % YES, CONSERVATIVE
            [Big_phi phi prob_cell MIP] = big_phi_comp_fb(M,x0_s,p,b_table,options,BRs,FRs);
        else % NO, PROGRESSIVE
            [Big_phi phi prob_cell MIP] = big_phi_comp_fb(M,x0_s,p,b_table,options);
        end
        MIP_M{M_i} = MIP;
    else
        [Big_phi phi prob_cell] = big_phi_comp(M,x0_s,p,b_table,options);
    end
    Big_phi_M(M_i) = Big_phi; % Big_phi for each subset
    phi_M{M_i} = phi; % Set of small_phis for each subset of each subset
    
    % I THINK THE POINTS IN CONCEPT SPACE ARE SOMEHOW ENCODED IN THE CELL
    % ARRAY prob_M
    prob_M{M_i,1} = prob_cell{1}; % prob_cell{2} & prob_cell{3} are repertoires in parts <--?
    prob_M{M_i,2} = prob_cell{2};
end


%% display
fprintf('\n')
fprintf('--------------------------------------------------------------\n')
fprintf('Big phi values in subset M\n')
for M_i = 1: 2^N-1
    fprintf('M=%s: Big_phi=%f\n',mod_mat2str(M_cell{M_i}),Big_phi_M(M_i));
end