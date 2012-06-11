function M2 = pick_rest(M,M1)

M2 = M;
for i=1: length(M1)
    M2(M2==M1(i)) = [];
end
