function prob = partial_prob_comp_time(T,x0_so,x0_in,x0_out,x1_b,source,p,b_table,op_fb)

% T: time interval over which phi is computed

if T == 1
    prob = partial_prob_comp2(x0_so,x0_in,x0_out,x1_b,source,p,b_table,op_fb);
elseif T == 2
    part = [x0_so x0_in];
    part = sort(part);
    N_p = length(part);
    if op_fb == 0 || op_fb == 2
        % forward
        N_t = length(x1_b);
       
        prob1 =  partial_prob_comp2(x0_so,x0_in,x0_out,part,source,p,b_table,op_fb); % (2^N_p, 1)
        
        prob2 = zeros(2^N_t,2^N_p);
        x0_so = part;
        x0_in = [];
        for j=1: 2^N_p
            x0 = b_table{j,N_p};
            source(part) = x0;
            prob2(:,j) = partial_prob_comp2(x0_so,x0_in,x0_out,x1_b,source,p,b_table,op_fb);
        end
        prob = prob2*prob1;
    else
        % backward
        N_t = length(x0_so);

        prob1 =  partial_prob_comp2(part,[],x0_out,x1_b,source,p,b_table,op_fb); % (2^N_p,1)
        
        prob2 = zeros(2^N_t,2^N_p);
        for j=1: 2^N_p
            x1 = b_table{j,N_p};
            source(part) = x1;
            prob2(:,j) = partial_prob_comp2(x0_so,x0_in,x0_out,part,source,p,b_table,op_fb);
        end
        prob = prob2*prob1;
    end
else
    
end


end


% function prob = partial_prob_comp_nn()
    