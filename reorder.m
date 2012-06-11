function x_re = reorder(M,x)
x_re = x;
for i=1: length(x)
    x_re(i) = find(M==x(i));
end
end