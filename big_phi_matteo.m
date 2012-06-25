function big_phi = big_phi_matteo(C, C_phi, steps_per_dimension)

% C = N+1 x Nc matrix of concepts as columns, N = dimension of simplex, Nc
%   is the number of concepts

% C_phi = 1xNc vector of phi values for each concept

% steps_per_dimension is how many points are seached in each dimension

% big_phi = BIG PHI!!

% get dimension of simplex
N = size(C,1) - 1;
% get number of concepts
Nc = size(C,2);

radius = 2/(N+3);

%% Build R

% initialize R, a matrix of coordinates of the simplex in the projected
% space where each column is the coordinates of one of the vertices of the
% simplex

R = zeros(N,N+1);
R(1,1) = 1;
R(1,2:end) = -1/N;

for j = 2:N
    
    R(j,j) = nthroot(sum(R(:,j).^N),N);
    R(j,j+1:end) = -1/N - R(1:j-1,j)'*R(1:j-1,j+1);

end

R(N,N+1) = -R(N,N);

%% Build the transform matrix T which maps from concept space to the space of R

r0_mat = repmat(R(:,1),1,N);
T = inv(R(:,2:end) - r0_mat);



%% Test all points in lattice bounded by [-1,1] for all dimensions with granularity passed in arguments

% numPoints = (2/grain + 1)^N;
% test_points = cell(numPoints,1);

ref_point = ones(N,1) * -1;
distance = zeros(N,1);

nSimplex = 0; big_phi_unnormalized = 0;
step_size = 2 / (steps_per_dimension - 1);

center_point = ones(N+1,1)/(N+1);
nPerSphere = 0;

for i = 0:(steps_per_dimension^N-1)
    
   for j = 1:N
       
       distance(j) = mod(i,steps_per_dimension);
       i = (i - distance(j)) / steps_per_dimension;
      

   end
   
%        disp('DISTANCE');
%        disp(distance);
       
       test_point = ref_point + distance * step_size;
%        disp('TEST POINT:');
%        disp(test_point)
%        
       r0 = zeros(N,1); r0(1) = 1;

       gammaTemp = T * (test_point - r0);
       gamma = zeros(N+1,1);
       gamma(1) = 1 - sum(gammaTemp);
       gamma(2:end) = gammaTemp;

       if (sum(gamma) == 1)
           nSimplex = nSimplex + 1;
           
           if(sum(abs(gamma - center_point)) < radius)
               nPerSphere = nPerSphere + 1;
           end

           phi_contribution = 0;
           for k = 1:size(C,2)

               if(sum(abs(gamma - C(:,k))) < radius && C_phi(k) > phi_contribution)
                   phi_contribution = C_phi(k);
               end

           end
           
           

           big_phi_unnormalized = big_phi_unnormalized + phi_contribution;
       
       end
       
end


% big_phi = big_phi_unnormalized * 5^N / (nSimplex * (3^N - N - 1));
big_phi = big_phi_unnormalized / nPerSphere;

disp(big_phi)


end



function [nSimplex big_phi_unnormalized] = test(test_point, nSimplex, big_phi_unnormalized, T, C, C_phi, N)

r0 = zeros(N,1); r0(1) = 1;

gammaTemp = T * (test_point - r0);
gamma = zeros(N+1,1);
gamma(1) = 1 - sum(gammaTemp);
gamma(2:end) = gammaTemp;

if (sum(gamma) ~= 1)
    return
end

nSimplex = nSimplex + 1;

phi_contribution = 0;
for i = 1:size(C,2)
   
    if(sum(abs(gamma - C(:,i))) < radius && C_phi(i) > phi_contribution)
        phi_contribution = C_phi(i);
    end
        
end
    
big_phi_unnormalized = big_phi_unnormalized + phi_contribution;



end



% an attempt at a recursive solution
% function [nSimplex big_phi_unnormalized] = test_points(current_point,grain,dim,N, nSimplex, big_phi_unnormalized,C, C_phi,T)
% 
% 
% if(dim < N || current_point(dim+1) ~= 1)
%     
%     move = zeros(N,1);
%     move(dim+1) = grain;
%     
%     next_point = current_point + move;
%     
%     [nSimplex big_phi_unnormalized] = test_points(next_point,grain,dim+1,N, nSimplex, big_phi_unnormalized,C,C_phi,T);
% end
% 
% disp(current_point);
% 
% % get point in simplex space
% r0 = zeros(N,1); r0(1) = 1;
% 
% gammaTemp = T * (current_point - r0);
% gamma = zeros(N+1,1);
% gamma(1) = 1 - sum(gammaTemp);
% gamma(2:end) = gammaTemp;
% 
% if (sum(gamma) ~= 1)
%     return
% end
% 
% nSimplex = nSimplex + 1;
% 
% phi_contribution = 0;
% for i = 1:size(C,2)
%    
%     if(sum(abs(gamma - C(:,i))) < radius && C_phi(i) > phi_contribution)
%         phi_contribution = C_phi(i);
%     end
%         
% end
%     
% big_phi_unnormalized = big_phi_unnormalized + phi_contribution;
% 
% end







    
    



