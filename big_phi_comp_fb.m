function [Big_phi phi prob_cell MIP M_IRR] = big_phi_comp_fb(M,x0_s,p,options)

%%  compute big phi for a subset, M
% M: a subset of the whole system (can be the whole system itself)
% x0_s: given data about the current state
% p: transition probability matrix in the whole system (p(x1|x0))

% THE FINAL TWO ARGS ARE OPTIONS, IF THEY ARE THERE THEN WE ARE DOING
% CONSERVATIVE, OTHERWISE PROGRESSIVE...

global grain
global J
% global BRs, global FRs
    
N = length(M);

op_fb = options(1);
op_phi = options(2);
op_figures = options(3);
op_single = options(4);
op_context = options(6);
op_min = options(9);
op_console = options(10);
op_big_phi = options(11);

op_whole = 1;

if op_figures == 0 || op_figures == 1
    disp_flag = 0;
elseif op_figures == 2 && N ~= log2(size(p,1)) % or we are not dealing with the whole complex
    disp_flag = 0;
else
    disp_flag = 1;
    if op_fb == 0
        fig_p = 200;
    else
        fig_p = 100;
    end
end

%% x0 data

% ???This is where we build numerators of purviews (power-set exclude empty
% set)
C_x0 = cell(2^N-1,1);
k = 1;
for i = 1:N % can this be done in one for-loop over k = 1:2^N-1 ?
    C = nchoosek(M,i); % create a matrix of combinations of M of size i
    N_C = size(C,1);
    % fprintf('i=%d N_c=%d\n',i,N_C);
    for j = 1:N_C % for all combos of size i
        x0 = C(j,:); % pick a combination
        C_x0{k} = x0;% store combo
        k = k + 1;
    end
end
M_p = C_x0; % power set of M

MIP = cell(2^N-1,1); % MIP in the past, the present, and the future <--?
phi_all_values = zeros(2^N-1,3); % small phis, for each purview, overall, backward, and forward

prob = cell(2^N-1,1); % transition repertoire
prob_prod = cell(2^N-1,1); % partitioned transition repertoire

M_IRR = cell(0,0);

%% computing small phis
EmptyCon = zeros(2^N-1,1);
parfor ci=1: 2^N-1  % loop over purview numerators
    x0 = C_x0{ci}; % given data of x0
    %Larissa smart purviews
    %if any element inside the numerator does not have inputs or outputs,
    %no need to calculate purview
    Nconnect = [sum(J(x0,:),2) sum(J(:,x0))'];
    EmptyCon(ci) = numel(Nconnect)-nnz(Nconnect);
    %EmptyCon(ci) =0; % Old version
    if EmptyCon(ci) == 0
%         if op_console
%             fprintf('C=%s\n',mod_mat2str(x0));
%         end
        if op_context == 0
            [phi_all_values(ci,:) prob{ci} prob_prod{ci} MIP{ci}] ...
            =  phi_comp_ex(options,M,x0,x0_s,p,M_p,J);
        else
            [phi_all_values(ci,:) prob{ci} prob_prod{ci} MIP{ci}] ...
                =  phi_comp_ex(options,M,x0,x0_s,p,M_p);
        end
    else
        phi_all_values(ci,:) = [0 0 0];
        prob{ci} = []; prob_prod{ci} = []; MIP{ci} = []; % should we change these to uniform, full sys... etc
    end    
end

prob_cell = cell(2,1);
prob_cell{1} = prob;
prob_cell{2} = prob_prod;


phi = phi_all_values(:,1);
phi_m = zeros(N,3); % cumulative sum

% PRETTY SURE THIS CAN JUST BE DONE WITH A SUM() CALL
for i_C=1: 2^N-1
    C = C_x0{i_C};
    i = length(C);
    phi_m(i,1) = phi_m(i,1) + phi(i_C);
    phi_m(i,2) = phi_m(i,2) + phi(i_C)/nchoosek(N,i);
end

% THIS SEEMS LIKE IT CAN BE DONE A SMARTER WAY
for i=1: N
    if i > 1
        phi_m(i,3) = phi_m(i-1,3) + phi_m(i,1);
    else
        phi_m(i,3) = phi_m(i,1);
    end
end


if (op_big_phi == 0)

    Big_phi = phi_m(end,3);


elseif (op_big_phi == 1)


    index_vec_IRR = find(phi ~= 0);
    N_IRR = length(index_vec_IRR);

    if(N_IRR~=0)

        concepts = zeros(2^N,N_IRR);
        concept_phis = zeros(1,N_IRR);


        j = 1;
        for i = 1:2^N-1

            if (phi(i) ~= 0)

                concepts(:,j) = prob{i}{1};
                concept_phis(j) = phi(i);
                j = j + 1;
            end

        end

        Big_phi = big_phi_volume(concepts,concept_phis,grain);

    else
        Big_phi = 0;
    end
elseif (op_big_phi == 2 || op_big_phi == 4)

    index_vec_IRR = find(phi ~= 0);
    N_IRR = length(index_vec_IRR);

    if(N_IRR~=0)

        concepts = zeros(2^N,N_IRR);
        concept_phis = zeros(1,N_IRR);


        j = 1;
        for i = 1:2^N-1

            if (phi(i) ~= 0)

                concepts(:,j) = prob{i}{1};
                concept_phis(j) = phi(i);
                j = j + 1;
            end

        end

        M_IRR = cell(N_IRR,1);

        for i=1: N_IRR
            j = index_vec_IRR(i);
            M_IRR{i} = C_x0{j};
        end

        if (op_big_phi == 2)
            Big_phi = big_phi_info(M_IRR,concepts,concept_phis);
        else
           Big_phi = NaN;
        end

    else
        Big_phi = 0;
    end
elseif (op_big_phi == 3)



    index_vec_IRR = find(phi ~= 0);
    N_IRR = length(index_vec_IRR);

    if(N_IRR > 1)

        concepts_past = zeros(2^N,N_IRR);
        concepts_future = zeros(2^N,N_IRR);
        concept_phis = zeros(1,N_IRR);


        j = 1;
        for i = 1:2^N-1

            if (phi(i) ~= 0)

                concepts_past(:,j) = prob{i}{1};
                concepts_future(:,j) = prob{i}{2};
                concept_phis(j) = phi(i);
                j = j + 1;

            end

        end

        if (N == size(p,2))
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

for i_C=1: 2^N-1
    C = C_x0{i_C};
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
                fprintf('\tPast: phi_b = %f\n',KLD(prob{i_C}{1},prob_prod{i_C}{1}));
%                 fprintf('Partition: %s\n',string_p{3});
%                 fprintf('Distribution (full):\n');
%                 disp(prob{i_C}{1}');
%                 fprintf('Distribution (partitioned):\n');
%                 disp(prob_prod{i_C}{1}');
                fprintf('\tFuture: phi_f = %f\n',KLD(prob{i_C}{2},prob_prod{i_C}{2}));
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

% prob_cell2 = cell(2^N-1,2);
% for i=1: 2^N-1
%     prob_cell2{i,1} = prob{i};
%     prob_cell2{i,2} = prob_prod{i};
% end


% if disp_flag == 1
%     if op_console
%         if size(p,2) == N
%             fprintf('\n')
%             fprintf('------------------------------------------------------------------------\n')
%             fprintf('Whole system\n')
%             fprintf('Core concepts: MIP: Small phi\n')
%         end
%     end
% %     plot_REP(prob_cell2,phi,MIP,100, M, op_context, op_min);
% end
% 
% if op_console
%     for i=1: N
%         fprintf('%d: phi_cum=%f phi_sum=%f phi_mean=%f\n',i,phi_m(i,3),phi_m(i,1),phi_m(i,2));
%     end
% 
%     fprintf('Big phi=%f\n',Big_phi);
% end

end
