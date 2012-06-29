function h = discrete_entropy(p)

p(p == 0) = [];

h = -sum(p.*log2(p));