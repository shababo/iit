function prob = partial_prob_comp_bfp(M,xp_1,xp_2,x0_1,x0_2,xf_1,xf_2,x0_s,p,b_table,op_whole,op_context)
% computing the partitioned transition repertoire

% M: subsystem
% xp,x0,xf: units inside the perspective at the past, the present and the
% future
% x0_s: given data of x0
% p: transition probability matrix (TPM)
% op_whole: 0: only targets 1: the whole system
% op_context: 0: conservative 1: progressive

if nargin < 11
    op_whole = 0;
end

N = size(p,2); % total number of the units in the whole system
N_M = length(M); % total number of the units in the subset M

xp = combine(xp_1,xp_2);
x0 = combine(x0_1,x0_2);
xf = combine(xf_1,xf_2);

Np = length(xp);
N0 = length(x0);
Nf = length(xf);

M_out = pick_rest(1:N,M); % units outside of the subsystem
x0_in = pick_rest(M,x0);
N0_in = length(x0_in);

%% partitioned backward
xp_in = pick_rest(M,xp);
xp_out = M_out;
xp1_out = pick_rest(1:N,xp_1);
xp2_out = pick_rest(1:N,xp_2);

xp_re = reorder(M,xp);
xp_1_re = reorder(xp,xp_1);
xp_2_re = reorder(xp,xp_2);

%% partitioned forward
xf_out = pick_rest(M,xf);
x0_out2 = combine(M_out,x0);

xf_re = reorder(M,xf);
xf_out_re = reorder(M,xf_out);
xf_1_re = reorder(xf,xf_1);
xf_2_re = reorder(xf,xf_2);

x01_so = pick_rest(M,x0_2);
x01_out = combine(x0_2,M_out);
x01_out2 = pick_rest(1:N,x0_1);
x02_so = pick_rest(M,x0_1);
x02_out = combine(x0_1,M_out);
x02_out2 = pick_rest(1:N,x0_2);

%% small system
if op_whole == 0
    p_b = zeros(2^Np,2^N0_in);
    p_f = zeros(2^Nf,2^N0_in);
    
    for i=1: 2^N0_in
        if N0_in ~= 0
            x0_s(x0_in) = b_table{i,N0_in};
        end
        % backward
        for j=1: 2^Np
            xp_s = b_table{j,Np};
            % outside of the perspective
            if op_context == 0
                prob_x0_in =  partial_prob_forward([],[],1:N,x0_in,[],x0_s(x0_in),p,b_table); % all out
            else
                prob_x0_in =  partial_prob_forward(xp,xp_in,xp_out,x0_in,xp_s,x0_s(x0_in),p,b_table);
            end
            p_b(j,i) = prob_x0_in;
            
            % partition 1
            if isempty(xp_1) == 0
                prob_x0_p1 =  partial_prob_forward(xp_1,[],xp1_out,x0_1,xp_s(xp_1_re),x0_s(x0_1),p,b_table);
                p_b(j,i) = p_b(j,i)*prob_x0_p1;
            end
            % partition2
            if isempty(xp_2) == 0
                prob_x0_p2 =  partial_prob_forward(xp_2,[],xp2_out,x0_2,xp_s(xp_2_re),x0_s(x0_2),p,b_table);
                p_b(j,i) = p_b(j,i)*prob_x0_p2;
            end
        end
        
        % forward
        for j=1: 2^Nf
            xf_s = b_table{j,Nf};
            
            p_f(j,i) = 1;
            % partition 1
            if isempty(xf_1) == 0
                if op_context == 0
                    prob_xf_p1 = partial_prob_forward(x0_1,[],x01_out2,xf_1,x0_s(x0_1),xf_s(xf_1_re),p,b_table);
                else
                    prob_xf_p1 = partial_prob_forward(x01_so,[],x01_out,xf_1,x0_s(x01_so),xf_s(xf_1_re),p,b_table);
                end
                p_f(j,i) = p_f(j,i)*prob_xf_p1;
            end
            
            % partition 2
            if isempty(xf_2) == 0
                if op_context == 0
                    prob_xf_p2 = partial_prob_forward(x0_2,[],x02_out2,xf_2,x0_s(x0_2),xf_s(xf_2_re),p,b_table);
                else
                    prob_xf_p2 = partial_prob_forward(x02_so,[],x02_out,xf_2,x0_s(x02_so),xf_s(xf_2_re),p,b_table);
                end
                p_f(j,i) = p_f(j,i)*prob_xf_p2;
            end
            
        end
        
        % fprintf('x0=%s p_b=%s p_f=%s\n',mat2str(x0_s),mat2str(p_b(:,i)),mat2str(p_f(:,i)));
    end
    
%% whole system
else
    
    p_b = zeros(2^N_M,2^N0_in);
    p_f = zeros(2^N_M,2^N0_in);
    
    for i=1: 2^N0_in
        if N0_in ~= 0
            x0_s(x0_in) = b_table{i,N0_in};
        end
        % backward
        for j=1: 2^N_M
            xp_whole = b_table{j,N_M};
            xp_s = xp_whole(xp_re);
            % outside of the perspective, prob_x0_in
            if op_context == 0
                p_b(j,i) =  partial_prob_forward([],[],1:N,x0_in,[],x0_s(x0_in),p,b_table);
            else
                p_b(j,i) =  partial_prob_forward(xp,xp_in,xp_out,x0_in,xp_s,x0_s(x0_in),p,b_table);
            end
            % partition 1, prob_x0_p1
            if isempty(xp_1) == 0
                p_b(j,i) =  p_b(j,i)*partial_prob_forward(xp_1,[],xp1_out,x0_1,xp_s(xp_1_re),x0_s(x0_1),p,b_table);
            end
            % partition2, prob_x0_p2
            if isempty(xp_2) == 0
                p_b(j,i) =  p_b(j,i)*partial_prob_forward(xp_2,[],xp2_out,x0_2,xp_s(xp_2_re),x0_s(x0_2),p,b_table);
            end
        end
        
        % forward
        for j=1: 2^N_M
            xf_whole = b_table{j,N_M};
            xf_s = xf_whole(xf_re);
            xf_s_out = xf_whole(xf_out_re);
            
            % outside of the perspective, prob_xf_out
            if op_context == 0
                p_f(j,i) = partial_prob_forward([],[],1:N,xf_out,[],xf_s_out,p,b_table); % all out
            else
                p_f(j,i) = partial_prob_forward(x0_in,[],x0_out2,xf_out,x0_s(x0_in),xf_s_out,p,b_table);
            end
            % partition 1, prob_xf_p1
            if isempty(xf_1) == 0
                if op_context == 0
                    p_f(j,i) = p_f(j,i)*partial_prob_forward(x0_1,[],x01_out2,xf_1,x0_s(x0_1),xf_s(xf_1_re),p,b_table);
                else
                    p_f(j,i) = p_f(j,i)*partial_prob_forward(x01_so,[],x01_out,xf_1,x0_s(x01_so),xf_s(xf_1_re),p,b_table);
                end
            end
            % partition 2, prob_xf_p2
            if isempty(xf_2) == 0
                if op_context == 0
                    p_f(j,i) = p_f(j,i)*partial_prob_forward(x0_2,[],x02_out2,xf_2,x0_s(x0_2),xf_s(xf_2_re),p,b_table);
                else
                    p_f(j,i) = p_f(j,i)*partial_prob_forward(x02_so,[],x02_out,xf_2,x0_s(x02_so),xf_s(xf_2_re),p,b_table);
                end
            end
            
        end
    end
    
end

prob = p_b*p_f';
Norm = sum(sum(prob));

if Norm ~= 0
    prob = prob/Norm;
end

end