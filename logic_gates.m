function output = logic_gates(input,logic_type)
% LOGIC_GATES the probability an element will turn on given the inputs
%
% OUTPUT = logic_gates(INPUT, LOGIC_TYPE)
%
% Given a binary input vector, input, and the logic type, logic_type, of
% the element of interest, this function will return the probability that
% the element of interest will be on


% AND
if logic_type == 1
    
    output = all(input);

% OR
elseif logic_type == 2

    output = any(input);

% XOR
elseif logic_type == 3

    output = sum(input) == 1;
    
% COPY    
elseif logic_type == 4

    output = input(1);

% NOT - TODO: check that element only has one input
elseif logic_type == 5
    
    output = ~input(1);

% NULL
elseif logic_type == 6

    output = .5;
end