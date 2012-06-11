function [BRs FRs] = comp_pers(x0_s,p,b_table,options) 

%%  compute BRs and FRs in every possible perspectives

N = size(p,2);

M_pb = cell(2^N,1); % based on binary counting
for i=1: 2^N
    temp = b_table{i,N};
    M_pb{i} = find(temp==1);
end

BRs = cell(2^N,2^N); % backward repertoire
FRs = cell(2^N,2^N); % forward repertoire
for i=1: 2^N
    x0 = M_pb{i}; % current
    x0_out = pick_rest(1:N,x0);
    x0_si = x0_s(x0);
    p_b_cell = cell(2^N-1,1);
    p_f_cell = cell(2^N-1,1);
    parfor j=1:2^N
        xp = M_pb{j}; % past
        xp_out = pick_rest(1:N,xp);
        p_b_cell{j} = partial_prob_cons(xp,[],xp_out,x0,x0_si,p,b_table,1); % all out except for xp
        
        xf = M_pb{j}; % future
        p_f_cell{j} = partial_prob_cons(x0,[],x0_out,xf,x0_si,p,b_table,0); % all out except for x0
    end
    
    % does this need to be here because they aren't distributed arrays?
    % ugh, also this for loop is TOTALLY unnecessary!
%     BRs{i,:} = p_b_cell{:};
%     FRs{i,:} = p_b_cell{:};
    for j=1: 2^N
        BRs{i,j} = p_b_cell{j};
        FRs{i,j} = p_f_cell{j};
    end
end