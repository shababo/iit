function big_phi(varargin)

%% check that there are either no arguments or all arguments
if (nargin ~= 0 && nargin ~= 3)
    fprintf('\nYou have not entered a valid number of arguments...\n\n');
    fprintf(['This function can be run with no arguments, in which case\n'...
             'you will be prompted for information about the system,\n'...
             'or you must run big_phi(J,logic_types,state),\n'...
             'where J is the connectivity matrix, logic_types is a vector\n'...
             'of logical types, and state is the current state of the system.\n']);
end

%% if no arguments, get information about system from user
if(nargin == 0)
    
    fprintf('\n');
    N = input('How many nodes in the system?  ');
    
    logic_types = zeros(N,1);
    
    fprintf(['\nPlease enter logic types for each element.\n'...
             '1-> AND, 2-> OR, 3-> XOR, 4 -> COPY, 5-> NOT, 6 -> NULL\n']);
    for i = 1:N
        valid = 0;
        while(~valid)
            type = input(['Logic type for element ' num2str(i) ':  ']);
            if (type > 0 && type < 7)
                valid = 1;
            else
                fprintf('Invalid logic type, please enter again...\n');
            end
        end
        logic_types(i) = type;     
    end
    
    fprintf(['\nPlease enter connectivity for each element.\n'...
             'Enter connectivity as a vector of efferents. For example,\n'...
             'if element 1 has directed edges to 2 and 3, enter: "[2 3]"\n']);
    
    J = zeros(N,N);
    for i = 1:N
        valid = 0;
        while(~valid)
            efferents = input(['Connectivity for element ' num2str(i) ':  ']);
            if (all(efferents > 0) && all(efferents <= N))
                valid = 1;
            else
                fprintf('Invalid vector, please enter again...\n');
            end
        end
        J(i,efferents) = 1;     
    end
    
    fprintf(['\nPlease enter the state of each element.\n'...
             'If an element is OFF enter 0 and enter 1 for ON.\n']);
         
    current_state = zeros(N,1);
    
    for i = 1:N
        valid = 0;
        while(~valid)
            state = input(['State for element ' num2str(i) ':  ']);
            if (state == 0 || state == 1)
                valid = 1;
            else
                fprintf('Invalid state, please enter either 0 or 1...\n');
            end
        end
        current_state(i) = state;     
    end
    
    disp(logic_types)
    disp(J)
    disp(current_state)

%%
else
    
end

