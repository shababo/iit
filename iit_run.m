function iit_run(tpm, in_J, current_state, in_noise, options)


N = size(tpm,2);

global grain, global noise, global BRs, global FRs, global J, global b_table

global output_data



grain = 50;
noise = in_noise;
BRs = cell(2^N,2^N); % backward repertoire
FRs = cell(2^N,2^N); % forward repertoire

J = in_J;

output_data.tpm = tpm;
output_data.J = J;
output_data.current_state = current_state;
output_data.noise = noise;
output_data.options = options;
output_data.num_nodes = N;

op_ave = options(18);

if op_ave == 0
    z_max = 1;
else
    z_max = 2^N;
end

% if nargin == 4 
%     J = ones(N);
% elseif nargin == 5
%     J = in_J;
% end



% binary table and states list
% from now on, all loops over subsets/bipartitions will use
% b_table for their ordering
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

output_data.states = states;

% parallel computing
op_parallel = options(19);
isOpen = matlabpool('size');
if  isOpen == 0 && op_parallel > 0
%     s = ['matlabpool ' int2str(op_parallel)];
%     eval(s);
    matlabpool;
end

op_complex = options(15);
op_fb = options(1);
op_context = options(6);

% init cell arrays for results
Big_phi_M_st = cell(2^N,1);
Big_phi_MIP_st = cell(2^N,1);
MIP_st = cell(2^N,1);
Complex_st = cell(2^N,1);
prob_M_st = cell(2^N,1);
phi_M_st = cell(2^N,1);
concept_MIP_M_st = cell(2^N,1);
complex_MIP_M_st = cell(2^N,1);


for z=1: z_max
    
    if op_ave == 0
        x1 = current_state;
    else
        x1 = states(:,z); 
    end
    
%     fprintf('x1=%s\n',mat2str(x1));
    
    % partial_prob_comp(partition, partition, state, prob_matrix, binary
    % table, op_fb
    check_prob = partial_prob_comp(1:N,1:N,x1,tpm,b_table,1); % last argument is op_fb = 1;
    state_check = sum(check_prob);
    if state_check == 0
        fprintf('This state cannot be realized!\n')
        Big_phi_M_st{z} = NaN;
        Big_phi_MIP_st{z} = NaN;
    else
        if op_complex == 0 % only consider whole system
            if op_fb == 2
                options(1) = 0; [Big_phi_f phi_f prob_cell_f] = big_phi_comp(1:N,x1,tpm,b_table,options);
                options(1) = 1; [Big_phi_b phi_b prob_cell_b] = big_phi_comp(1:N,x1,tpm,b_table,options);
                Big_phi = Big_phi_f + Big_phi_b;
            elseif op_fb == 3 % THIS IS THE ONLY ONE WE DO NOW? BOTH FORWARD AND BACKWARD SIMULTANEOUSLY
                M = 1:N;
                if op_context == 0
%                     [BRs FRs] = comp_pers(x1,tpm,b_table,options);
                    [Big_phi phi prob_cell MIPs M_IRR] = big_phi_comp_fb(M,x1,tpm,b_table,options);
                    % irreducible points
                    [IRR_REP IRR_phi IRR_MIP M_IRR] = IRR_points(prob_cell,phi,MIPs,M, 0,op_fb);
                    fprintf('\n')
                    fprintf('---------------------------------------------------------------------\n\n')
                    fprintf('Big_phi = %f\n', Big_phi);
                    fprintf('Sum of small_phis = %f\n',sum(phi));
                    fprintf('\nCore Concepts For Complex (Purview, MIP(past & future), Small phi):\n\n');
                    plot_REP(Big_phi, IRR_REP,IRR_phi,IRR_MIP, 1, M, options)
                else
                    [Big_phi phi prob_cell MIP prob_cell2] = big_phi_comp_fb(M,x1,tpm,b_table,options);
                end
            else
                [Big_phi phi prob_cell] = big_phi_comp(1:N,x1,tpm,b_table,options);
            end
            Big_phi_st(z) = Big_phi;
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % THE CURRENT SETTINGS TAKE US HERE    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else % find the complex
            
            [MIP Complex Big_phi_M Big_phi_MIP_M prob_M phi_M concept_MIP_M complex_MIP_M M_cell] ...
                = big_phi_complex(x1,tpm,options)
            
            if op_fb == 2
%                 % subindex b means backward and f means forward
%                 IRR_phi_b = IRR_phi{1};
%                 IRR_phi_f = IRR_phi{2};
%                 IRR_REP_b = IRR_REP{1};
%                 IRR_REP_f = IRR_REP{2};
%                 M_IRR_b = M_IRR{1};
%                 M_IRR_f = M_IRR{2};
%                 % state dependent big phi and big phi MIP
%                 Big_phi_st(z) = sum(IRR_phi_b)+sum(IRR_phi_f);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % THE CURRENT SETTINGS TAKE US HERE    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                Big_phi_M_st{z} = Big_phi_M;
                Big_phi_MIP_st{z} = Big_phi_MIP_M;
                MIP_st{z} = MIP;
                Complex_st{z} = Complex;
                prob_M_st{z} = prob_M;
                phi_M_st{z} = phi_M;
                concept_MIP_M_st{z} = concept_MIP_M;
                complex_MIP_M_st{z} = complex_MIP_M;
                
            end
        end
    end
    
    % pause;
end

% load up output_data struct
output_data.Big_phi_M = Big_phi_M_st;
output_data.Big_phi_MIP = Big_phi_MIP_st;
output_data.MIP = MIP_st;
output_data.Complex = Complex_st;
output_data.concepts_M = prob_M_st;
output_data.small_phi_M = phi_M_st;
output_data.concept_MIP_M = concept_MIP_M_st;
output_data.complex_MIP_M = complex_MIP_M_st;
output_data.M_cell = M_cell;

% if op_ave == 1
%     if op_fb == 0
%         Big_phi_ave = sum(Big_phi_st)/2^N;
%     else
%         Big_phi_ave = sum(p_x1 .* Big_phi_st); %weighted ave/expected value
%     end
% %     fprintf('Big_phi_ave=%f\n',Big_phi_ave);
% end

op_close = 0;
isOpen = matlabpool('size');
if isOpen > 0 && op_close == 1
    matlabpool close;
end

fprintf('Loading GUI... \n');

save('output_data_sample.mat','output_data');

iit_explorer(output_data)