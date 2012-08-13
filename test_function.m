function test_function

state_size_vec = [2 4 5 2];

% check vectorized function
tic
for i=0:prod(state_size_vec)-1
    dec2multibase_array(i,state_size_vec);
end
toc
% check for-loop function
tic
for i = 0:prod(state_size_vec)-1
    dec2multibase_matlab(i,state_size_vec);
end
toc
% check mex
tic
for i=0:prod(state_size_vec)-1
    dec2multibase(i,state_size_vec)'
end
toc
% check first two functions match
for i=0:prod(state_size_vec)-1
    
    if(any(dec2multibase_array(i,state_size_vec)' ~= dec2multibase_matlab(i,state_size_vec)))
        disp(['error on ' num2str(i)])
        disp(dec2multibase_array(i,state_size_vec)')
        disp(dec2multibase(i,state_size_vec))
    end
    
end
% check last two functions match
for i=0:prod(state_size_vec)-1
    
    if(any(dec2multibase_matlab(i,state_size_vec) ~= dec2multibase(i,state_size_vec)))
        disp(['error on ' num2str(i)])
        disp(dec2multibase_array(i,state_size_vec)')
        disp(dec2multibase(i,state_size_vec))
    end
    
end

