function [small_phi concept_points concept_MIP MIP] = small_phi(tpm, present_set, cur_state, op_direction, op_find_purviews, op_phi_function, varargin)

% small_phi = small_phi(tpm, present_set, cur_state, op_direction, op_find_purviews, op_phi_function, [past_set, future_set, J])
% 
% Input:
%
% tpm                   state X node transition probability matrix
% present_set           set of elements in the numerator of the purview
% cur_state             state of all nodes in the system in the present
% op_direction          compute backwards (0), forwards (1), or both (2)
% op_find_purview       find the purview (0) by choosing max(phi_mip) or
%                       compute small_phi for the given purview (1) (see last 
%                       two arguments)
% op_phi_function       use KLD (0) or EMD (1)
% past_set              for op_find_purview == 1, the set of nodes to look
%                       at in the past
% future_set            for op_find_purview == 1, the set of nodes to look
%                       at in the future
% J                     optional connectivity matrix
%
% Output:
% small_phi             the small_phi value of the purview

global BRs, global FRs


% Fill in unset optional values.
switch nargin
    case 6
        past_set = [];
        future_set = [];
        J = [];
        op_find_purviews = 0;
    case 7
        if op_direction == 0
            past_set = varargin{1};
        elseif op_direction == 1
            future_set = varargin{1};
        end
    case 8
        if op_direction == 0
            past_set = varargin{1};
            J = varargin{2};
        elseif op_direction == 1
            future_set = varargin{1};
            J = varargin{2};
        elseif op_direction == 2
            past_set = varargin{1};
            future_set = varargin{2};
            J = [];
        end
    case 9
        past_set = varargin{1};
        future_set = varargin{2};
        J = varargin{3};
end

N = size(tpm,2);
system_set = 1:N;

BRs = cell(2^N,2^N); % backward repertoire
FRs = cell(2^N,2^N); % forward repertoire

% binary table
b_table = cell(2^N,N);
states = zeros(N,2^N);
for i=1: N
    for j=1: 2^i
        b_table{j,i} = trans2(j-1,i);
        if i== N
            states(:,j) = trans2(j-1,i);
        end
    end
end

options = zeros(1,16);
options(6) = 0; % op_context
options(7) = 0; % op_whole
options(9) = 1; % op_min
options(14) = 1; % op_normalize
options(16) = op_phi_function; % KLD or EMD

if (op_find_purviews == 0)
    
    % build the power set
    power_set = cell(2^N-1,1);
    k = 1;
    for i = 1:N % can this be done in one for-loop over k = 1:2^N-1 ?
        C = nchoosek(system_set,i); % create a matrix of combinations of M of size i
        N_C = size(C,1);

        for j = 1:N_C % for all combos of size i
            x0 = C(j,:); % pick a combination
            power_set{k} = x0;% store combo
            k = k + 1;
        end
    end

    phi_MIP = zeros(2^N-1,2);
    prob_cand = cell(2^N-1,1);
    prob_prod_MIP_cand = cell(2^N-1,1);
    MIP_cand = cell(2^N-1,1);
    
    for i=1: 2^N-1
        
        %Larissa smart purviews: Only test those connections that actually exist
        other_set = power_set{i};
        
        if ~isempty(J) && nnz(sum(J(present_set,other_set),1) == 0) > 0 % some x is not input of x0 (numerator) --> no phiBR
            if nnz(sum(J(other_set,present_set),2) == 0) == 0 % but x is output
                [phi_MIP(i,:) prob_cand{i} prob_prod_MIP_cand{i} MIP_cand{i}] ...
                    = phi_comp_bORf(options,present_set,other_set,tpm,2,b_table,cur_state); 
            else
                prob_cand{i} = cell(2,1);
                prob_prod_MIP_cand{i} = cell(2,1);
                MIP_cand{i} = cell(2,2,2);
            end
        else
            if ~isempty(J) && nnz(sum(J(other_set,present_set),2) == 0) > 0 % x is not output, but x is input
                [phi_MIP(i,:) prob_cand{i} prob_prod_MIP_cand{i} MIP_cand{i}] ...
                    = phi_comp_bORf(options,present_set,other_set,tpm,1,b_table,cur_state); 
            else % x is both or J is empty
                [phi_MIP(i,:) prob_cand{i} prob_prod_MIP_cand{i} MIP_cand{i}] ...
                    = phi_comp_bf(options,system,present_set,other_set,other_set,cur_state,tpm,b_table); 
            end 
        end    
    end
    
    % exlusion principle
    max_phi_MIP_bf = zeros(2,1); % backward and forward phi
    MIP = cell(2,2,2);
    concept_points = cell(2,1);
    concept_MIP = cell(2,1);
    for bf = 1:2
        [max_phi_MIP_bf(bf) j_max] = max_ex(phi_MIP(:,bf),power_set);
        MIP(:,:,bf) = MIP_cand{j_max}(:,:,bf);
        concept_points{bf} = prob_cand{j_max}{bf};
        concept_MIP{bf} = prob_prod_MIP_cand{j_max}{bf};
        if bf == 1
            past_set = power_set{j_max};
        else
            future_set = power_set{j_max};
        end
    end
    
    if (op_direction == 1 || op_direction == 2)
       small_phi = max_phi_MIP_bf(op_direction);

    else
       small_phi = min(max_phi_MIP_bf(1),max_phi_MIP_bf(2));
    end
    
    
else
    
    
    [phi_MIP concept_points concept_MIP MIP] ...
    = phi_comp_bf(options,system,present_set,past_set,future_set,cur_state,tpm,b_table);
    
    if (op_direction == 1 || op_direction == 2)
       small_phi = phi_MIP(op_direction);

    else
       small_phi = min(phi_MIP(1),phi_MIP(2));
    end
    
end

%% imposing maxent on units outside of perspectives

for i = 1:2
    if i == 1
        x = past_set;
    else
        x = future_set;
    end
    if length(x) ~= N
        concept_points{i} = expand_prob(concept_points{i},system,x);
        concept_MIP{i} = expand_prob(concept_MIP{i},system,x);
    end
end
