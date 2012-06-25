% b_table = a cell array of binary values (2^N x N)

% BRs = a cell array of

% f = firing rate of each element (N x 1) NB: assumes that all states are
%       equally represented during normal firing patterns

% J = connectivy matrix (N x N)

% logic_type = type of logic gate for each element (N x 1); 
%       KEY: logic type: 1-> AND, 2-> OR, 3-> XOR, 4 -> COPY, 5-> NOT,
%       6 -> NULL

% M_cell = a cell array of subsets

% N = number of elementary units in the system

% Na = number of afferents... not sure why you need this... couldn't it be
%       different for each element?

% p = the TPM... why not call it TPM??!?! in any case, this is currently a
%       (2^N x N) matrix which means that it tells you the state of each
%       element at t+1 given the state at t

% p_x0 = given a state, what are the
%       values of each element at the next state (2^N x N)

% p_x1 = the probability of each state

% p1 = transition matrix - given the current state, what is the prob of each next state
%       (2^N x 2^N)

% states = a stored table of binary values of length N (N x 2^N)
 
% z_max = number of states to compute over (z_max = 1 when you onl look at
%       the current state, z_max = 2^N when you compute (and average) over
%       all possible states