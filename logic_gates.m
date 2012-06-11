function output = logic_gates(i_vec,logic_type)

if logic_type== 1
    % AND gate
    if i_vec(1) == 1 && i_vec(2)==1
        output = 1;
    else
        output = 0;
    end
elseif logic_type == 2
    % OR gate
    if i_vec(1) == 0 && i_vec(2)==0
        output = 0;
    else
        output = 1;
    end
elseif logic_type == 3
    % XOR gate
    if i_vec(1) == i_vec(2)
        output = 0;
    else
        output = 1;
    end
elseif logic_type == 4
    % COPY gate
    output = i_vec(1);
elseif logic_type == 5
    % NOT gate
    output = 1 - i_vec(1);
elseif logic_type == 6
    % NULL gate
    output = 1/2;
end