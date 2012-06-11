op = 1; % 0:Difference of entropy 1:KL-divergence

N = 4;
Na = 3;
op_parallel = 0;
op_close = 0;

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

%% x0 bipartion
N_b = 0; % number of bipartition
for i=1: N-1
    N_b = N_b + nchoosek(N,i);
end
x0_b1_o = cell(N_b,1);
x0_b2_o = cell(N_b,1);
i_b = 1;
for i=1: N-1
    C_b = nchoosek(1:N,i);
    N_C = size(C_b,1);
    for j=1: N_C
        x0_b1_o{i_b} = C_b(j,:);
        temp = 1: 1: N;
        temp(C_b(j,:)) = [];
        x0_b2_o{i_b} = temp;
        i_b = i_b + 1;
    end
end

%% x1 data
C_x1 = cell(2^N-1,1);
k = 1;
for i=1: N
    C = nchoosek(1:N,i);
    N_C = size(C,1);
    H = zeros(N_C,1);
    phi = zeros(N_C,1);
    % fprintf('i=%d N_c=%d\n',i,N_C);
    for j=1: N_C
        C_j = C(j,:);
        if i>1
            [B1 B2 N1_b] = bipartition(C_j);
        end
        C_x1{k} = C_j;
        k = k + 1;
    end
end

isOpen = matlabpool('size');
if  isOpen == 0 && op_parallel > 0
    s = ['matlabpool ' int2str(op_parallel)];
    eval(s);
end

%% Big phi
Big_phi_x1 = zeros(2^N,1);
p_x1 = zeros(2^N,1);

for z=2^N: 2^N % 2^N
    x1 = trans2(z-1,N);
    % x1 = [0; 0; 0; 1; 1; 1];
    display(x1)
    for k=1: 2^N
        p_x1(z) = p_x1(z) + sig_prob(x1,1:N,k,p);
    end
    p_x1(z) = p_x1(z)/2^N;
    % display(p_x1(z))
    % pause;
    
    MIP0 = cell(2^N-1,1);
    MIP1 = cell(2^N-1,1);
    
    H = zeros(2^N-1,1);
    phi = zeros(2^N-1,1);
    parfor i_C=1: 2^N-1
        fprintf('i_C=%d\n',i_C);
        C_j = C_x1{i_C}; % given data of x1
        i = length(C_j);
        prob = ones(2^N,1);
        for k=1: 2^N
            prob(k) = sig_prob(x1,C_j,k,p);
        end
        if sum(prob) ~= 0
            prob = prob/sum(prob);
            prob(prob==0) = 1;
            
            H(i_C) = -sum(prob.*log2(prob)); % conditional entropy in the whole system
            
            display(C_j);
            % display(prob);
            fprintf('H in whole=%f\n',H(i_C));
            % pause;
            
            % computing the conditional entropy in the parts considering only bipartition
            
            if i > 1
                [phi(i_C) MIP0{i_C} MIP1{i_C}] =  phi_comp(x0_b1_o,x0_b2_o,C_j,x1,N_b,H(i_C),p,op);
            else
                phi(i_C) = N - H(i_C);
            end
            
            
        end
        
        % pause;
        
    end
    
    phi_m = zeros(N,3);
    
    for i_C=1: 2^N-1
        C = C_x1{i_C}; 
        i = length(C);
        fprintf('C%d=[',i);
        for k=1: i
            fprintf('%d ',C(k))
        end
        fprintf(']: ');
        
        if i > 1
            x0_p1 = MIP0{i_C};
            x1_p1 = MIP1{i_C};
            
            fprintf('[')
            for k=1: length(x0_p1)
                fprintf('%d ',x0_p1(k));
            end
            fprintf(']-[');
            for k=1: length(x1_p1)
                fprintf('%d ',x1_p1(k));
            end
            fprintf('] ');
        end
        
        if abs(phi(i_C)) < 10^-8
            phi(i_C) = 0;
        end
        fprintf('phi%d=%f\n',i,phi(i_C));
        
        phi_m(i,1) = phi_m(i,1) + phi(i_C);
        phi_m(i,2) = phi_m(i,2) + phi(i_C)/nchoosek(N,i);
        
    end
    
    for i=1: N
        if i > 1
            phi_m(i,3) = phi_m(i-1,3) + phi_m(i,1);
        else
            phi_m(i,3) = phi_m(i,1);
        end
        fprintf('%d: phi_cum=%f phi_sum=%f phi_mean=%f\n',i,phi_m(i,3),phi_m(i,1),phi_m(i,2));
    end
    
    Big_phi = phi_m(end,3);
    fprintf('Big_phi=%f\n',Big_phi);
    
    figure(1)
    subplot(1,3,1),bar(1:N,phi_m(:,1))
    subplot(1,3,2),bar(1:N,phi_m(:,2))
    subplot(1,3,3),bar(1:N,phi_m(:,3))
    
    figure(2)
    bar(phi)
    
    Big_phi_x1(z) = Big_phi;
    
end

Big_phi_ave = sum(p_x1.*Big_phi_x1);

 isOpen = matlabpool('size');
if isOpen > 0 && op_close == 1
matlabpool close;
end