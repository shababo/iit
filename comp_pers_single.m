function rep = comp_pers_single(current,other,x0_s,p,bf_option)

%  compute BRs and FRs for a single perspective but given some fixed
%  current state

global b_table nodes


N = size(p,2);

% x0 = find(b_table{current,N} == 1);
x0 = current;

% x0 = M_pb{current}; % current subset
x0_out = pick_rest(1:N,x0); % complement
x0_si = x0_s(x0); % get the current state of the current subset

if (bf_option == 1) % backwards
%     xp = M_pb{past}; % past
%     xp = find(b_table{other,N} == 1);
    xp = other;
    xp_out = pick_rest(1:N,xp);
    rep = partial_prob_cons(xp,[],xp_out,x0,x0_si,p,b_table,1); % all out except for xp
    
    if ~isempty(xp)
        matrix_size = ones(1,N);
        matrix_size(xp) = 2;
        rep = reshape(rep,matrix_size);
    end
    
    
end

if (bf_option == 2) % forwards
%     xf = M_pb{future}; % future
%     xf = find(b_table{other,N} == 1);
    xf = other;
    rep = partial_prob_cons(x0,[],x0_out,xf,x0_si,p,b_table,0); % all out except for x0
    
    if ~isempty(xf)
        matrix_size = ones(1,N);
        matrix_size(xf) = 2;
        rep = reshape(rep,matrix_size);
    end
end



    

end