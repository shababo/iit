function [Big_phi_MIP MIP Big_phi_M IRR_phi IRR_REP states] = big_phi_MIP(x1,p,b_table,options)

%% compute Big-phi and find MIP
% x1 : given data about the current neural state
% p: transition matrix in the whole system
% op: the way how small phi is computed (0: difference of entropy, 1: KLD)
% op_disp: (1: display the results, 0: not)

N = log2(size(p,1)); % number of elements in the whole system

%% subset
M_cell = cell(2^N-1,1);
M_cell2 = cell(2^N-1,1);

% k = 1;
% for i=1: N
%     C = nchoosek(1:N,i);
%     N_C = size(C,1);
%     for j=1: N_C
%         C_j = C(j,:);
%         M_cell{k} = C_j;
%         k = k + 1;
%     end
% end

for i=1: 2^N-1
    x = trans2(i,N);
    C = [];
    for j=1: N
        if x(j) == 1
            C = [C j];
        end
    end
    M_cell{i} = C;
end

%% compute big phi in every subset
Big_phi_M = zeros(2^N-1,1);
phi_M = cell(2^N-1,1);
prob_M = cell(2^N-1,1);
for M_i = 1: 2^N-1
    M = M_cell{M_i};
    
    fprintf('M=%s\n',mod_mat2str(M));
    [Big_phi phi phi_m prob] = big_phi_comp(M,x1,p,b_table,options);
    Big_phi_M(M_i) = Big_phi;
    phi_M{M_i} = phi;
    prob_M{M_i} = prob;
end

%% Find Big-phi MIP in the whole system
% [Big_phi_MIP MIP] = MIP_search(1:N,N,Big_phi_M);

%% Find complex
op_complex = 1;
if op_complex == 1
    Big_phi_MIP_M = zeros(2^N-1,1);
    MIP_M = cell(2^N-1,1);
    
    % Big_phi_MIP_M(end) = Big_phi_MIP;
    % MIP_M{end} = MIP;
    % index_vec = sort_index(N);
    for M_i = 1: 2^N-1
        % M_i_d = index_vec(M_i+1)-1;
        % M = M_cell{M_i_d};
        M = M_cell{M_i};
        % M_cell2{M_i} = M;
        % fprintf('M_i=%d M_i_d=%d M=%s\n',M_i,M_i_d,mod_mat2str(M));
        if length(M) > 1
            [Big_phi_MIP_M(M_i) MIP_M{M_i}] = MIP_search(M,N,Big_phi_M);
        end
    end
    
    [Big_phi_MIP_comp M_i_max] = max(Big_phi_MIP_M);
    % Complex = M_cell2{M_i_max};
    Complex = M_cell{M_i_max};
    MIP_comp = MIP_M{M_i_max};
    M2 = pick_rest(Complex,MIP_comp);
    fprintf('Complex: %s MIP=%s-%s Big_phi_MIP=%f\n', ...
        mat2str(Complex),mod_mat2str(MIP_comp),mod_mat2str(M2),Big_phi_MIP_comp);
    
    IRR_REP = prob_M{M_i_max};
    IRR_phi = phi_M{M_i_max};
    IRR_REP(:,IRR_phi==0) = [];
    IRR_phi(IRR_phi==0) = [];
    IRR_phi = IRR_phi';
    [index_vec states] = sort_index(length(Complex));
    IRR_REP = IRR_REP(index_vec,:);
    
    Big_phi_MIP = Big_phi_MIP_comp;
    MIP = MIP_comp;
    
    op_disp = options(3);
    if op_disp ~= 0
        fig_max = 8;
        N_IRR = length(IRR_phi);
        for i_C = N_IRR: 1
            figure(3)
            subplot(fig_max,1,N_ir-fig_max*fig_pi),bar(prob_all(index_vec))
            axis([-Inf Inf 0 y_max]);
            sC = make_title(C);
            sPhi = [sC,': \phi=',num2str(phi(i_C),4)];
            title(sPhi)
            if i_C == 1
                sy = ['\Phi=',num2str(Big_phi,4)];
                xlabel(sy)
            end
        end
    end
    
%   options(end) = options(end)*3;
%   [Big_phi_comp IRR_phi IRR_phi_m IRR_REP_cell] = big_phi_comp(Complex,x1,p,b_table,options);
%     IRR_REP_cell(IRR_phi==0) = [];
%     IRR_phi(IRR_phi==0) = [];
%     N_IRR = length(IRR_phi);
%     IRR_REP = zeros(2^length(Complex),N_IRR);
%     for i=1: N_IRR
%         IRR_REP(:,i) = IRR_REP_cell{i};
%     end
%     IRR_phi = IRR_phi';
%     [index_vec states] = sort_index(length(Complex));
%     IRR_REP = IRR_REP(index_vec,:);
end