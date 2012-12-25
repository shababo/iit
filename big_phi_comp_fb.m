function [Big_phi phi_all_values prob_cell MIP M_IRR network] = big_phi_comp_fb(subsystem,whole_sys_state,network)

%%  compute big phi for a subset, subsystem
% subsystem: a subset of the whole system (can be the whole system itself)
% whole_sys_state: given data about the current state
% p: transition probability matrix in the whole system (p(x1|numerator))

% THE FINAL TWO ARGS ARE OPTIONS, IF THEY ARE THERE THEN WE ARE DOING
% CONSERVATIVE, OTHERWISE PROGRESSIVE...

% global BRs, global FRs
num_nodes_subsys = length(subsystem);
num_states_subsys = prod([network.nodes(subsystem).num_states]);

op_single = network.options(4);
op_context = network.options(6);
op_min = network.options(9);
op_console = network.options(10);
op_big_phi = network.options(11);


%% numerator data

% ???This is where we build subsets_subsys of purviews (power-set exclude empty
% set)
subsets_subsys = cell(num_states_subsys - 1, 1);

% we can do this better for sure - TODO Larissa: Do it as in big_phi_all
% the subsets!!
k = 1;
for subset_size = 1:num_nodes_subsys % can this be done in one for-loop over k = 1:num_states_subsys-1 ?
    C = nchoosek(subsystem,subset_size); % create a matrix of combinations of subsystem of size i
    N_C = size(C,1);
    % fprintf('i=%d N_c=%d\n',i,N_C);
    for j = 1:N_C % for all combos of size i
        numerator = C(j,:); % pick a combination
        subsets_subsys{k} = numerator;% store combo
        k = k + 1;
    end
end

MIP = cell(num_states_subsys-1,1); % MIP in the past, the present, and the future <--?
phi_all_values = zeros(num_states_subsys-1,3); % small phis (for each purview) overall, backward, and forward

prob = cell(num_states_subsys-1,1); % transition repertoire
prob_prod = cell(num_states_subsys-1,1); % partitioned transition repertoire

M_IRR = cell(0,0);

%% computing small phis
EmptyCon = zeros(num_states_subsys-1,1);
for ci=1: num_states_subsys-1  % loop over purview subsets_subsys
    numerator = subsets_subsys{ci}; % given data of numerator
    %Larissa smart purviews
    %if any element inside the numerator does not have inputs or outputs,
    %no need to calculate purview
    Nconnect = [sum(network.connect_mat(numerator,:),2) sum(network.connect_mat(:,numerator))'];
    EmptyCon(ci) = numel(Nconnect)-nnz(Nconnect);
    %EmptyCon(ci) =0; % Old version
    if EmptyCon(ci) == 0
%         if op_console
%             fprintf('C=%s\n',mod_mat2str(numerator));
%         end
        [phi_all_values(ci,:) prob{ci} prob_prod{ci} MIP{ci} network] ...
            =  phi_comp_ex(subsystem,numerator,whole_sys_state,subsets_subsys,network);

    else
        phi_all_values(ci,:) = [0 0 0];
        uniform_dist = ones(1,num_states_subsys)/num_states_subsys;
        prob{ci} = {uniform_dist, uniform_dist}; 
        prob_prod{ci} = {uniform_dist, uniform_dist}; 
        MIP{ci} = cell(2,2,2); % should we change these to uniform, full sys... etc
    end    
end

prob_cell = cell(2,1);
prob_cell{1} = prob;
prob_cell{2} = prob_prod;


phi = phi_all_values(:,1);
phi_m = zeros(num_nodes_subsys,3); % cumulative sum

index_vec_IRR = find(phi ~= 0);
N_IRR = length(index_vec_IRR);

if(N_IRR~=0)

    concepts = zeros(num_states_subsys,N_IRR);
    concept_phis = zeros(1,N_IRR);


    j = 1;
    for i = 1:num_states_subsys-1

        if (phi(i) ~= 0)

            concepts(:,j) = prob{i}{1};
            concept_phis(j) = phi(i);
            j = j + 1;
        end

    end

    M_IRR = cell(N_IRR,1);

    for i=1: N_IRR
        j = index_vec_IRR(i);
        M_IRR{i} = subsets_subsys{j};
    end

end 
    


% PRETTY SURE THIS CAN JUST BE DONE WITH A SUM() CALL
for i_C=1: num_states_subsys-1
    C = subsets_subsys{i_C};
    i = length(C);
    phi_m(i,1) = phi_m(i,1) + phi(i_C);
    phi_m(i,2) = phi_m(i,2) + phi(i_C)/nchoosek(num_nodes_subsys,i);
end

% THIS SEEMS LIKE IT CAN BE DONE A SMARTER WAY
for i=1: num_nodes_subsys
    if i > 1
        phi_m(i,3) = phi_m(i-1,3) + phi_m(i,1);
    else
        phi_m(i,3) = phi_m(i,1);
    end
end


if (op_big_phi == 0 || op_big_phi == 6 || op_big_phi == 7) %Larissa: 6 is now L1 distance, was not assigned before anyways

    Big_phi = phi_m(end,3);

% Larissa: All that follows seems to be a mess. Not sure what is working at
% all
elseif (op_big_phi == 1)


%     index_vec_IRR = find(phi ~= 0);
%     N_IRR = length(index_vec_IRR);

    if(N_IRR~=0)

%         concepts = zeros(num_states_subsys,N_IRR);
%         concept_phis = zeros(1,N_IRR);
% 
% 
%         j = 1;
%         for i = 1:num_states_subsys-1
% 
%             if (phi(i) ~= 0)
% 
%                 concepts(:,j) = prob{i}{1};
%                 concept_phis(j) = phi(i);
%                 j = j + 1;
%             end
% 
%         end

        Big_phi = big_phi_volume(concepts,concept_phis,grain);

    else
        Big_phi = 0;
    end
elseif (op_big_phi == 2 ||op_big_phi == 4 || op_big_phi == 5)

    index_vec_IRR = find(phi ~= 0);
    N_IRR = length(index_vec_IRR);

    if(N_IRR~=0)

        concepts = zeros(num_states_subsys,N_IRR);
        concept_phis = zeros(1,N_IRR);


        j = 1;
        for i = 1:num_states_subsys-1

            if (phi(i) ~= 0)

                concepts(:,j) = prob{i}{1};
                concept_phis(j) = phi(i);
                j = j + 1;
            end

        end

        M_IRR = cell(N_IRR,1);

        for i=1: N_IRR
            j = index_vec_IRR(i);
            M_IRR{i} = subsets_subsys{j};
        end

        if (op_big_phi == 2)    %Larissa ?? So for 4 and 5 it is always NaN anyways? Why calculate something before?
            Big_phi = big_phi_info(M_IRR,concepts,concept_phis);
        else
            Big_phi = NaN;
        end

    elseif (op_big_phi == 2)
        Big_phi = 0;
    else
        Big_phi = NaN;
    end
elseif (op_big_phi == 3)



    index_vec_IRR = find(phi ~= 0);
    N_IRR = length(index_vec_IRR);

    if(N_IRR > 1)

        concepts_past = zeros(num_states_subsys,N_IRR);
        concepts_future = zeros(num_states_subsys,N_IRR);
        concept_phis = zeros(1,N_IRR);


        j = 1;
        for i = 1:num_states_subsys-1

            if (phi(i) ~= 0)

                concepts_past(:,j) = prob{i}{1};
                concepts_future(:,j) = prob{i}{2};
                concept_phis(j) = phi(i);
                j = j + 1;

            end

        end

        if (num_nodes_subsys == size(p,2))
            display = 1;
        else
            display = 0;
        end

%         disp(big_phi_spacing(concepts_past,concept_phis,0));
%         disp(big_phi_spacing(concepts_future,concept_phis,0));


        Big_phi = big_phi_spacing(concepts_past,concept_phis,0) + big_phi_spacing(concepts_future,concept_phis,0);

    elseif (N_IRR == 1)
        Big_phi = phi(1);
    else
        Big_phi = 0;
    end
else
    err = MException('Options:UnSetValue', ...
        'Option for method of computing Big Phi is incorrect');
    throw(err);
end



%% display

for i_C=1: num_states_subsys-1
    C = subsets_subsys{i_C};
    i = length(C);
    
    if abs(phi(i_C)) < 10^-8
        phi(i_C) = 0;
    end
        
    if op_console
        % Larissa smart purviews
        if EmptyCon(i_C) > 0
%             fprintf('C=%s: nodes lack input or output \n',mod_mat2str(C))
        else
            if i > 1 || op_single == 1        
                [string_p string] = make_title_fb(MIP{i_C},op_context,op_min);
                fprintf('C = %s: Core Concept: %s\n',mod_mat2str(C),string{3});
                fprintf('MIP = %s\n', string_p{3});
                fprintf('\tPast: phi_b = %f\n',phi_all_values(i_C,2));
%                 fprintf('Partition: %s\n',string_p{3});
%                 fprintf('Distribution (full):\n');
%                 disp(prob{i_C}{1}');
%                 fprintf('Distribution (partitioned):\n');
%                 disp(prob_prod{i_C}{1}');
                fprintf('\tFuture: phi_f = %f\n',phi_all_values(i_C,3));
%                 fprintf('Distribution (full):\n');
%                 disp(prob{i_C}{2}');
%                 fprintf('Distribution (partitioned):\n');
%                 disp(prob_prod{i_C}{2}');
            end

            fprintf('\tphi_mip = %f\n\n',phi(i_C));
        end
    % prob{i_C}
    % prob_prod{i_C}
    end
end

% prob_cell2 = cell(num_states_subsys-1,2);
% for i=1: num_states_subsys-1
%     prob_cell2{i,1} = prob{i};
%     prob_cell2{i,2} = prob_prod{i};
% end


% if disp_flag == 1
%     if op_console
%         if size(p,2) == num_nodes_subsys
%             fprintf('\n')
%             fprintf('------------------------------------------------------------------------\n')
%             fprintf('Whole system\n')
%             fprintf('Core concepts: MIP: Small phi\n')
%         end
%     end
% %     plot_REP(prob_cell2,phi,MIP,100, subsystem, op_context, op_min);
% end
% 
% if op_console
%     for i=1: num_nodes_subsys
%         fprintf('%d: phi_cum=%f phi_sum=%f phi_mean=%f\n',i,phi_m(i,3),phi_m(i,1),phi_m(i,2));
%     end
% 
%     fprintf('Big phi=%f\n',Big_phi);
% end

end
