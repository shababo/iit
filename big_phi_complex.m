function [MIP Complex Big_phi_M Big_phi_MIP_M prob_M phi_M concept_MIP_M complex_MIP_M M_cell] = big_phi_complex(x0,p,options)
% function [Big_phi_MIP MIP Complex M_i_max Big_phi_M Big_phi_MIP_M prob_M phi_M concept_MIP_M complex_MIP_M M_cell] = big_phi_complex(x0,p,options)

% [Big_phi_MIP MIP Big_phi_M IRR_phi IRR_REP IRR_MIP M_IRR prob_M phi_M MIP_M] = big_phi_complex(x0,p,b_table,options)
% 
% INPUTS:
% x0            current state of the whole system as a binary Nx1 column vector
% p             transition probability matrix
% b_table       a cell array of binary values stored as vectors
% options       algorithm options
%
% OUTPUTS:
% Big_phi_MIP   Big phi MIP in the complex    
% MIP           MIP in the complex
% Big_phi_M     Big phi in every subsets
% IRR_phi       small phi values of irreducible points
% IRR_REP       actual repertoires of irreducible points
% IRR_MIP
% M_IRR
% prob_M
% phi_M
% MIP_M

global b_table
global output_data

%% options
op_fb = options(1); % 0: forward repertoire, 1: backward repertoire, 2: both, 3: simultaneous
op_phi = options(2); % two versions of small phi 0:Difference of entropy, 1:KL-divergence 
op_figures = options(3);  % 0: No figures, 1: only complex 2: complex and whole system, 3: all figures
op_console = options(10);
op_big_phi = options(11);
op_sum = options(12);

op_context = options(6);
op_min = options(9);

%%

N = size(p,2); % number of elements in the whole system

if op_fb == 2
    %forward computation
    options(1) = 0;
    [Big_phi_M_f phi_M_f prob_M_f] = big_phi_all(x0,p,b_table,options);
    % backward computation
    options(1) = 1;
    [Big_phi_M_b phi_M_b prob_M_b M_cell] = big_phi_all(x0,p,b_table,options);
    % sum up the forward and backward phi
    Big_phi_M = Big_phi_M_f + Big_phi_M_b;
    % complex search
    [Big_phi_MIP MIP Complex M_i_max] = complex_search(Big_phi_M,M_cell,N);
    % irreducible points
    [IRR_REP_f IRR_phi_f M_IRR_f] = IRR_points(prob_M_f,phi_M_f,Complex,M_i_max,op_fb);
    [IRR_REP_b IRR_phi_b M_IRR_b] = IRR_points(prob_M_b,phi_M_b,Complex,M_i_max,op_fb);
    
    % store both the forward and backward results
    IRR_REP = cell(2,1);
    IRR_phi = cell(2,1);
    M_IRR = cell(2,1);
    IRR_REP{1} = IRR_REP_f;
    IRR_REP{2} = IRR_REP_b;
    IRR_phi{1} = IRR_phi_f;
    IRR_phi{2} = IRR_phi_b;
    M_IRR{1} = M_IRR_f;
    M_IRR{2} = M_IRR_b;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE CURRENT SETTINGS TAKE US HERE    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif op_fb == 0 || op_fb == 1 || op_fb == 3
    % 0: forward or 1: backward computation or 3: simultaneous forward and
    % backward computation
    [Big_phi_M phi_M prob_M M_cell concept_MIP_M M_IRR_M] = big_phi_all(x0,p,options)
    % complex search
    [Big_phi_MIP MIP Complex M_i_max  Big_phi_MIP_M complex_MIP_M] = complex_search(Big_phi_M,M_cell, M_IRR_M, N,prob_M,phi_M,options);
    % irreducible points - do this outside of here...
%     [IRR_REP IRR_phi IRR_MIP M_IRR] = IRR_points(prob_M,phi_M,MIP_M,Complex, M_i_max,op_fb);
%     fprintf('Sum of small_phis = %f\n',sum(IRR_phi));
%     fprintf('\nCore Concepts For Complex (Purview, MIP(past & future), Small phi):\n\n');
%     plot_REP(Big_phi_M(M_i_max), IRR_REP,IRR_phi,IRR_MIP, 1, Complex, options)
end

%% plot irreducible points in the complex
% % op_figures = 1;
% if op_figures ~= 0
%     if op_fb == 2
%         if length(IRR_phi{1}) > 8 || length(IRR_phi{2}) > 8
%             fig_co = 2;
%         else
%             fig_co = 1;
%         end
%     else
%         fig_co = 1;
%     end
%     if op_fb == 2
%         for i=1: 2
%             plot_IRR(IRR_REP{i},IRR_phi{i},M_IRR{i}, op_fb+i-1, fig_max*fig_co,fig_co)
%         end
%     elseif op_fb == 0 || op_fb == 1
%         plot_IRR(IRR_REP,IRR_phi,M_IRR, op_fb, fig_max*fig_co,fig_co,IRR_MIP)
% %     elseif op_fb == 3
% %         fprintf('\nCore Concepts For Complex (Purview, MIP(past & future), Small phi):\n\n');
% % %         fprintf('Purview || MIP || Small phi\n');
% %         plot_REP(Big_phi_M(M_i_max), IRR_REP,IRR_phi,IRR_MIP, 1, Complex, op_context, op_min)
%         
%         
%     end
% end
% 
% 
% 
% 
% %% Plot irreducible points
% function [] = plot_IRR(IRR_REP,IRR_phi,M_IRR, op, fig_max, fig_co,IRR_MIP)
% N_IRR = length(IRR_phi);
% y_max = 0;
% 
% if op ~= 4
%     for i_C=1: N_IRR
%         y_max = max(y_max,max(IRR_REP(:,i_C)));
%     end
%     y_max = y_max + 0.02;
% end
% 
% reverse_vec = N_IRR: -1: 1;
% for i = 1: N_IRR
%     j = reverse_vec(i);
%     
%     fig_pi = floor((i-1)/fig_max);
%     fig_i = i - fig_pi*fig_max;
%     figure(1+fig_pi);
%     if op == 0 ||op == 1 || op == 4
%         if fig_co == 1
%             convert_vec = 1: fig_max;
%         else
%             convert_vec = [2*(1: fig_max/2)-1 2*(1: fig_max/2)];
%         end
%         if op== 4
%             sq = sqrt(fig_max);
%             prob = IRR_REP{j,1};
%             prob_prod = IRR_REP{j,2};
%             subplot(sq,sq,convert_vec(fig_i)),imagesc(prob);
%             colormap('gray')
%             caxis([0 1])
%             figure(11+fig_pi);
%             subplot(sq,sq,convert_vec(fig_i)),imagesc(prob_prod);
%             colormap('gray')
%             caxis([0 1])
%         else
%             subplot(fig_max/fig_co,fig_co,convert_vec(fig_i)),bar(IRR_REP(:,j))
%         end
%     elseif op == 2
%         if fig_co == 1
%             convert_vec = 2*(1: fig_max);
%         else
%             convert_vec = [4*(1: fig_max/2)-1 4*(1: fig_max/2)];
%         end
%         subplot(fig_max/fig_co,2*fig_co,convert_vec(fig_i)),bar(IRR_REP(:,j))
%     elseif op == 3
%         if fig_co == 1
%             convert_vec = 2*(1: fig_max) - 1;
%         else
%             convert_vec = [4*(1: fig_max/2)-3 4*(1: fig_max/2)-2];
%         end
%         subplot(fig_max/fig_co,2*fig_co,convert_vec(fig_i)),bar(IRR_REP(:,j))
%     end
%     if op~= 4
%         axis([-Inf Inf 0 y_max]);
%     end
%     
%     C = M_IRR{j};
%     if op == 4
%         [string_p string] = make_title_fb(IRR_MIP{j});
%         sC = string;
%         figure(11+fig_pi);
%         title(string_p)
%     else
%         op_s = mod(op,2);
%         sC = make_title(C,op_s);
%     end
%     
%     figure(1+fig_pi);
%     sPhi = [sC,': \phi=',num2str(IRR_phi(j),3)];
%     title(sPhi)
%     if i== 1
%         fprintf('Irreducible points\n');
%     end
%     fprintf('%s\n',sPhi);
%     
%     if j == 1
%         Big_phi_comp = sum(IRR_phi);
%         sy = ['\Phi=',num2str(Big_phi_comp,3)];
%         xlabel(sy)
%     end
%     
% end
