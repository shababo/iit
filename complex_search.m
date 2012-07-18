function [Big_phi_MIP MIP Complex M_i_max Big_phi_MIP_M MIP_M] = complex_search(Big_phi_M,M_cell,M_IRR_M,N,prob_M, phi_M,options)
%% Find complex

Big_phi_MIP_M = zeros(2^N-1,1);
MIP_M = cell(2^N-1,1);

fprintf('\n')
fprintf('\nBig_phi_MIP in subset M:\n\n')
for M_i = 1: 2^N-1
    M = M_cell{M_i};
    if length(M) > 1
        
        [Big_phi_MIP_M(M_i) MIP_M{M_i}] = MIP_search(M,N,Big_phi_M, M_IRR_M, prob_M, phi_M,options);
        
    end
end

[Big_phi_MIP M_i_max] = max_complex(Big_phi_MIP_M,M_cell);
Complex = M_cell{M_i_max};

disp('**********************')
disp(M_i_max == convi(Complex));
disp('**********************')

MIP = cell(2,1);
MIP{1} = MIP_M{M_i_max};
Big_phi = Big_phi_M(M_i_max);

MIP{2} = pick_rest(Complex,MIP{1});

fprintf('\n')
fprintf('---------------------------------------------------------------------\n\n')
 fprintf('Complex = %s\nBig_phi = %f\nMIP = %s-%s\nBig_phi_MIP = %f\n', ...
     mat2str(Complex),Big_phi, mod_mat2str(MIP{1}),mod_mat2str(MIP{2}),Big_phi_MIP);
 
 
end

function [x_max i_max] = max_complex(x,M_cell)

N = length(x);
epsilon = 10^-10;
x_max = x(1);
i_max = 1;
for i=2: N
    dif = x_max - x(i);
    if abs(dif) < epsilon
        if length(M_cell{i-1}) < length(M_cell{i})
            x_max = x(i);
            i_max = i;
        end
    elseif dif < 0
        x_max = x(i);
        i_max = i;
    else
    end
end

end