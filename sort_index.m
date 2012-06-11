function [index_vec states] = sort_index(N)

index_vec = zeros(2^N,1);
states = zeros(2^N,N);

index_vec(1) = 1;
i_b = 2;
for i=1: N
    C_b = nchoosek(1:N,i);
    N_C = size(C_b,1);
    for j=1: N_C
        C = C_b(j,:);
        x = binarize(C,N);
        x_i = trans10(x);
        index_vec(i_b) = x_i;
        states(i_b,:) = x;
        i_b = i_b + 1;
    end
end

function b = binarize(C,N)

N_C = length(C);
b = zeros(N,1);

for i=1: N_C
    b(C(i)) = 1;
end