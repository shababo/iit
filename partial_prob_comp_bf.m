function prob = partial_prob_comp_bf(M,x0,xp,xf,x0_s,p,b_table,op_whole,op_context)

% computing the transition repertoire, which is the conditional probability distribution of the past (xp) and the future
% state (xf) given the current state x0, p(xp,xf|x0) 

% M: subsystem
% xp,x0,xf: units inside the perspective at the past, the present and the
% future
% x0_s: given data of x0
% p: transition probability matrix (TPM)
% op_whole: 0: only targets 1: the whole system
% op_context: 0: conservative 1: progressive

if nargin < 8
    op_whole = 0;
end

N = size(p,2); % total number of the units in the whole system
N_M = length(M); % total number of the units in the subsystem

Np = length(xp);
Nf = length(xf);
N0 = length(x0);

M_out = pick_rest(1:N,M); % units outside of the subsystem
x0_in = pick_rest(M,x0);
x0_c = pick_rest(1:N,x0);
N0_in = length(x0_in);

xp_in = pick_rest(M,xp);
xp_out_in = M_out;
xp_out = pick_rest(1:N,xp);

xf_out = pick_rest(M,xf);
x0_out2 = combine(M_out,x0);

xp_re = reorder(M,xp);
xf_re = reorder(M,xf);
xf_out_re = reorder(M,xf_out);

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
        
        if op_context == 0
            prob_x0_in =  partial_prob_forward([],[],1:N,x0_in,[],x0_s(x0_in),p,b_table); % all out
        else
            prob_x0_in =  partial_prob_forward(xp,xp_in,xp_out_in,x0_in,xp_s,x0_s(x0_in),p,b_table);
        end
        prob_x0 =  partial_prob_forward(xp,[],xp_out,x0,xp_s,x0_s(x0),p,b_table);
        
        p_b(j,i) = prob_x0*prob_x0_in;
    end
    
    % forward
    for j=1: 2^Nf
        xf_s = b_table{j,Nf};
        x0_so = M;
        x0_out = M_out;
        if op_context == 0
            p_f(j,i) =  partial_prob_forward(x0,[],x0_c,xf,x0_s(x0),xf_s,p,b_table); % all out except for x0
        else
            p_f(j,i) = partial_prob_forward(x0_so,[],x0_out,xf,x0_s(x0_so),xf_s,p,b_table);
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
            
            if op_context == 0
                prob_x0_in =  partial_prob_forward([],[],1:N,x0_in,[],x0_s(x0_in),p,b_table); % all out
            else
                prob_x0_in =  partial_prob_forward(xp,xp_in,xp_out_in,x0_in,xp_s,x0_s(x0_in),p,b_table);
            end
            prob_x0 =  partial_prob_forward(xp,[],xp_out,x0,xp_s,x0_s(x0),p,b_table); % all out except for xp
            
            p_b(j,i) = prob_x0*prob_x0_in;
        end
        
        % forward
        for j=1: 2^N_M
            xf_whole = b_table{j,N_M};
            xf_s = xf_whole(xf_re);
            xf_s_out = xf_whole(xf_out_re);
            x0_so = M;
            x0_out = M_out;
            
            if op_context == 0
                prob_xf_in =  partial_prob_forward(x0,[],x0_c,xf,x0_s(x0),xf_s,p,b_table); % all out except for x0
                prob_xf_out = partial_prob_forward([],[],1:N,xf_out,[],xf_s_out,p,b_table); % all out
            else
                prob_xf_in =  partial_prob_forward(x0_so,[],x0_out,xf,x0_s(x0_so),xf_s,p,b_table);
                prob_xf_out = partial_prob_forward(x0_in,[],x0_out2,xf_out,x0_s(x0_in),xf_s_out,p,b_table);
            end
            
            p_f(j,i) = prob_xf_in*prob_xf_out;
        end
        
        % fprintf('x0=%s p_b=%s p_f=%s\n',mat2str(x0_s),mat2str(p_b(:,i)),mat2str(p_f(:,i)));
    end
    
    
    
end

prob = p_b*p_f';
Norm = sum(sum(prob));


if Norm ~= 0
    prob = prob/Norm;
end