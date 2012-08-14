function test_function

% state_size_vec = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];
% 
% % check vectorized function
% tic
% for i=0:prod(state_size_vec)-1
%     dec2multibase_array(i,state_size_vec);
% end
% toc
% % check for-loop function
% tic
% for i = 0:prod(state_size_vec)-1
%     dec2multibase_matlab(i,state_size_vec);
% end
% toc
% % check mex
% tic
% for i=0:prod(state_size_vec)-1
%     dec2multibase(i,state_size_vec);
% end
% toc
% 
% % compare to trans2.mex
% tic
% for i=0:prod(state_size_vec)-1
%     trans2(i,state_size_vec,15);
% end
% toc
% % check first two functions match
% % for i=0:prod(state_size_vec)-1
% %     
% %     if(any(dec2multibase_array(i,state_size_vec)' ~= dec2multibase_matlab(i,state_size_vec)))
% %         disp(['error on ' num2str(i)])
% %         disp(dec2multibase_array(i,state_size_vec)')
% %         disp(dec2multibase(i,state_size_vec))
% %     end
% %     
% % end
% % % check last two functions match
% % for i=0:prod(state_size_vec)-1
% %     
% %     if(any(trans2(i,10) ~= dec2multibase(i,state_size_vec)))
% %         disp(['error on ' num2str(i)])
% %         disp(dec2multibase_array(i,state_size_vec)')
% %         disp(dec2multibase(i,state_size_vec))
% %     end
% %     
% % end


set = [1 2 3 4 5 6 7 8 9];
tic
ismember(1,set)
toc
tic
any(1 == set)
toc