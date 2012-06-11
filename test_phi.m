M = [1 2 3]; % subset
x0 = [3]; % current
xp = [1 2]; % past
xf = [2 3]; % future
[phi_MIP prob prob_prod_MIP MIP] = phi_comp_bf(options,M,x0,xp,xf,x0_s,p,b_table);
fprintf('phi_MIP=%f\n',phi_MIP);

figure(1)
subplot(1,2,1),imagesc(prob)
subplot(1,2,2),imagesc(prob_prod_MIP)
colormap('gray')
