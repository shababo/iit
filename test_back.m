op = 1; % 0:Difference of entropy 1:KL-divergence

N = 4;

J = zeros(N,N);
phi_cell = cell(N,2);

for i=1: N
    for j=1: N
        if i ~= j
            J(i,j) = 1;
        else
            % J(i,j) = 1;
        end
    end
end

% J = zeros(N,N);
% for i=1: N/2
%     i1 = 2*i-1;
%     i2 = 2*i;
%     J(i1:i2,i1:i2) = [0 1; 1 0];
% end

Na = 3;
J = zeros(N,N);
for i=1: N
    x = 1: N;
    x(i) = [];
    y = randsample(N-1,Na);
    z = x(y);
    J(i,z) = (N-1)/Na;
end

% J = J_3;
% J = J_2;

% J = J_fix;

a = 1.5;
T = 0.01;
% x1 = ones(N,1);

p = zeros(N,2^N);

for i=1: 2^N
    x0 = trans2(i-1,N);
    p(:,i) = (1+tanh((J*x0-a)/T))/2; % probability of turning on given the data x0
end

f = zeros(N,1);
for i=1: N
    f(i) = sum(p(i,:))/2^N;
end
avef = sum(f)/N;
display(J);
display(p);
fprintf('avef=%f\n',avef);
pause;

%% bipartion
N_b = 0; % number of bipartition
for i=1: N-1
    N_b = N_b + nchoosek(N,i);
end
x0_b1 = cell(N_b,1);
x0_b2 = cell(N_b,1);
i_b = 1;
for i=1: N-1
    C_b = nchoosek(1:N,i);
    N_C = size(C_b,1);
    for j=1: N_C
        x0_b1{i_b} = C_b(j,:);
        temp = 1: 1: N;
        temp(C_b(j,:)) = [];
        x0_b2{i_b} = temp;
        i_b = i_b + 1;
    end
end
%%

Big_phi_x1 = zeros(2^N,1);
p_x1 = zeros(2^N,1);

for z=2^N: 2^N % 2^N
    x1 = trans2(z-1,N);
    % x1 = [0; 0; 1; 1];
    display(x1)
    for k=1: 2^N
        p_x1(z) = p_x1(z) + sig_prob(x1,1:N,k,p);
    end
    p_x1(z) = p_x1(z)/2^N;
    % display(p_x1(z))
    % pause;
    
    MIP = cell(2^N-1-N,2);
    MIP_i = 1;
    for i=1: N
        C = nchoosek(1:N,i);
        N_C = size(C,1);
        H = zeros(N_C,1);
        phi = zeros(N_C,1);
        fprintf('i=%d N_c=%d\n',i,N_C);
        % pause;
        for j=1: N_C
            Norm = 0;
            C_j = C(j,:); % given data of x1
            prob = ones(2^N,1);
            for k=1: 2^N
                %             for l=1: i
                %                 if x1(C_j(l)) == 1
                %                     prob(k) = prob(k)*p(C_j(l),k);
                %                 else
                %                     prob(k) = prob(k)*(1-p(C_j(l),k));
                %                 end
                %             end
                %             prob(k)
                prob(k) = sig_prob(x1,C_j,k,p);
                % pause;
            end
            if sum(prob) ~= 0
                prob = prob/sum(prob);
                prob(prob==0) = 1;
                
                H(j) = -sum(prob.*log2(prob)); % conditional entropy in the whole system
                
                display(C_j);
                % display(prob);
                fprintf('H in whole=%f\n',H(j));
                % pause;
                
                
                % computing the conditional entropy in the parts
                % only consider bipartition
                
                if i > 1
                    [x1_b1 x1_b2 N1_b] = bipartition(C_j); % partition of x1
                    phi_cand = zeros(N_b*N1_b,2);
                    H_vec = zeros(N_b*N1_b,2);
                    i_H = 1;
                    for i_b = 1: N_b
                        x0_1 = x0_b1{i_b};
                        x0_2 = x0_b2{i_b};
                        
                        for j_b=1: N1_b
                            x1_1 = x1_b1{j_b};
                            x1_2 = x1_b2{j_b};
                            
                            if op == 0
                                [H1 prob1] = partial_ent(x0_1,x1_1,x1,p); % conditional entropy H(x0_1|x1_1)
                                
                                %                     display(x0_1);
                                %                     display(x1_1);
                                %                     display(prob1);
                                %                     fprintf('H1=%f\n',H1);
                                %
                                [H2 prob2] = partial_ent(x0_2,x1_2,x1,p); % conditional entropy H(x0_2|x1_2)
                                
                                %                     display(x0_2);
                                %                     display(x1_2);
                                %                     display(prob2);
                                %                     fprintf('H2=%f\n',H2);
                                
                                %                     if i== N
                                %                         pause;
                                %                     end
                                % sum
                                
                            else
                                [H1 prob1] = partial_ent(x0_1,x1_1,x1,p,C_j); % conditional entropy H(x0_1|x1_1)
                                [H2 prob2] = partial_ent(x0_2,x1_2,x1,p,C_j); % conditional entropy H(x0_2|x1_2)
                            end
                            H_vec(i_H,1) = H1;
                            H_vec(i_H,2) = H2;
                            phi_cand(i_H,1) = (H1 + H2) - H(j);
                            phi_cand(i_H,2) = phi_cand(i_H,1)/min(length(x0_1),length(x0_2));
                            i_H = i_H + 1;
                        end
                        
                    end
                    
                    [min_norm_phi i_phi_min] = min(phi_cand(:,2));
                    
                    
                    j_b = rem(i_phi_min,N1_b);
                    if j_b == 0
                        j_b = N1_b;
                    end
                    i_b = (i_phi_min-j_b)/N1_b+1;
                    display(x0_b1{i_b})
                    display(x1_b1{j_b})
                    display(x0_b2{i_b})
                    display(x1_b2{j_b})
                    phi(j)  =  phi_cand(i_phi_min,1); % difference of entropy gives phi
                    
                    MIP{MIP_i,1} = x0_b1{i_b};
                    MIP{MIP_i,2} = x1_b1{j_b};
                    
                    
                else
                    phi(j) = N - H(j);
                end
                
                if i > 1
                    fprintf('phi=%f H=%f minH_sum=%f H1=%f H2=%f\n',phi(j),H(j),phi(j)+H(j),H_vec(i_phi_min,1),H_vec(i_phi_min,2));
                else
                    fprintf('phi=%f H=%f\n',phi(j),H(j));
                end
                
            end
            
            if i>1
                MIP_i = MIP_i + 1;
            end
            
            % pause;
            
        end
        
        phi_cell{i,1} =phi;
        phi_cell{i,2} = C;
        
        
    end
    
    Big_phi = 0;
    phi_m = zeros(N,3);
    phi_all = zeros(2^N-1,1);
    
    l = 1;
    MIP_i = 1;
    for i=1: N
        phi = phi_cell{i,1};
        C = phi_cell{i,2};
        N_C = size(C,1);
        phi_all(l:l+N_C-1) = phi;
        for j=1: N_C
            fprintf('C%d=[',i);
            for k=1: i
                fprintf('%d ',C(j,k))
            end
            fprintf(']: ');
            
            if i > 1
                x0_p1 = MIP{MIP_i,1};
                x1_p1 = MIP{MIP_i,2};
                
                fprintf('[')
                for k=1: length(x0_p1)
                    fprintf('%d ',x0_p1(k));
                end
                fprintf(']-[');
                for k=1: length(x1_p1)
                    fprintf('%d ',x1_p1(k));
                end
                fprintf('] ');
                MIP_i = MIP_i + 1;
            end
            fprintf('phi%d=%f\n',i,phi(j));
        end
        phi_m(i,1) = sum(phi);
        phi_m(i,2) = mean(phi);
        if i > 1
            phi_m(i,3) = phi_m(i,3) + phi_m(i-1,3) + sum(phi);
        else
            phi_m(i,3) = phi_m(i,3) + sum(phi);
        end
        fprintf('%d: phi_cum=%f phi_sum=%f phi_mean=%f\n',i,phi_m(i,3),phi_m(i,1),phi_m(i,2));
        Big_phi = Big_phi + phi_m(i,1);
        
        l = l + N_C;
    end
    
    fprintf('Big_phi=%f\n',Big_phi);
    
    for i=1: N
        if abs(phi_m(i,1)) < 10^(-8)
            phi_m(i,1) = 0;
            phi_m(i,2) = 0;
        end
    end
    figure(1)
    subplot(1,3,1),bar(1:N,phi_m(:,1))
    subplot(1,3,2),bar(1:N,phi_m(:,2))
    subplot(1,3,3),bar(1:N,phi_m(:,3))
    
    figure(2)
    bar(phi_all)
    
    Big_phi_x1(z) = Big_phi;
    
end

Big_phi_ave = sum(p_x1.*Big_phi_x1);