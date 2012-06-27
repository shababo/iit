function [] = plot_REP(Big_phi, REP_cell,phi_vec,MIP_cell, fig_st, M, op_context, op_min)

% close all;
save bug_check;
% pause;

if nargin < 6
    op_context = 0; % default: conservative
end
if nargin < 7
    op_min = 1;
end

N = length(phi_vec);

figure(20+fig_st);
plot_phi(phi_vec,MIP_cell, M)
title('\phi  For Each Purview','FontSize',20);
xlabel('Purview','FontSize',16)
ylabel('\phi','FontSize',20)

if op_context == 0
    if N > 8
        r = 16;
        c = 2;
        fig_max = 16;
    else
        r = N;
        c = 2;
        fig_max = 8;
    end
    
    pw = 300; % panel width including margin in pixels
    ph = 50; % panel height including margin in pixels
    mb = 20; % bottom margin
    mt = 20; % top margin
    mh = 25; % margin height
    mw = 10; % margin width
    % panel_size = [0, 0, pw, ph]'; % panel size including margin
    panel_size_w = [50, 0, pw-mw*2, ph-mh]'; % panel size without margin
    fig_size = [0, 0,  pw*c + 100, ph*r + mb+mt + 50]'; % size of figure window
    pos_fig = [100,100,0,0]' + fig_size; % position of figure
    
    pos_vec = zeros(4,r,c);
    for i=1: r
        for j=1: c
            pos_vec(:,i,j) = [(j-1)*pw + 2*(j-1)+mw, (r-i)*ph + 10, 0, 0]' + panel_size_w + [0 mb 0 0]'; % pixels
            pos_vec([1 3],i,j) = pos_vec([1 3],i,j)/pos_fig(3); % normalization
            pos_vec([2 4],i,j) = pos_vec([2 4],i,j)/pos_fig(4); % normalization
        end
    end
    
%     Big_phi = sum(phi_vec);
    sy = ['Phi=',num2str(Big_phi)];
    
    for i=1 : N
        i_rev = N-i+1;
        fig_pi = floor((i_rev-1)/fig_max);
        fig_i = i_rev - fig_pi*fig_max;
        pos_BR = pos_vec(:,fig_i,1);
        pos_FR = pos_vec(:,fig_i,2);
        
        prob = REP_cell{i,1};
        prob_prod = REP_cell{i,2};
        
        BR_w = prob{1};
        FR_w = prob{2};
        BR_p = prob_prod{1};
        FR_p = prob_prod{2};
        phi_b = KLD(BR_w,BR_p);
        phi_f = KLD(FR_w,FR_p);
        [string_p string] = make_title_fb(MIP_cell{i},op_context,op_min);
        
        s_title = cell(2,1);
        s_title_p = cell(2,1);
%         s_title{1}= [string{1},': \phi_b=',num2str(phi_b)];
%         s_title{2}= [string{2},': \phi_f=',num2str(phi_f)];
%         s_title_p{1} = [string_p{1},': \phi_b=',num2str(phi_b)];
%         s_title_p{2} = [string_p{2},': \phi_f=',num2str(phi_f)];
        
        s_title{1}= {['FULL:    ', string{1}, ',  \phi_b=',num2str(phi_b)],['PARTITIONED:    ', string_p{1}]};
        s_title{2}= {['FULL:    ', string{2}, ',  \phi_f=',num2str(phi_f)],['PARTITIONED:    ', string_p{2}]};
        
        figure(fig_st+fig_pi)
        set(gcf,'Position',pos_fig)
        if mod(fig_i,min(r,N)) == 0
            labelON = 1;
        else
            labelON = 0;
        end
        BR = ones(length(BR_w),2);
        BR(:,1) = BR_w; BR(:,2) = BR_p;
        FR = ones(length(FR_w),2);
        FR(:,1) = FR_w; FR(:,2) = FR_p;
        plot_BRFR(BR,FR,pos_BR',pos_FR',s_title,labelON)
%         plot_BRFR(BR_w,FR_w,pos_BR',pos_FR',s_title,labelON)
        
        if i == N
            
            uicontrol('Style', 'text',...
            'String', {'Core Concept Distributions For Complex',sy},... %replace something with the text you want
            'Units','normalized',...
            'FontSize',16,...
            'BackgroundColor','w',...
            'Position', [0.33 0.9 0.33 0.08]);
        
            uicontrol('Style', 'text',...
            'String', 'PAST <-- CURRENT',... %replace something with the text you want
            'Units','normalized',...
            'FontSize',13,...
            'BackgroundColor','w',...
            'Position', [0.15 0.79 0.28 0.05]);
        
            uicontrol('Style', 'text',...
            'String', 'CURRENT --> FUTURE',... %replace something with the text you want
            'Units','normalized',...
            'FontSize',13,...
            'BackgroundColor','w',...
            'Position', [0.59 0.79 0.28 0.05]);
        
            my_legend = legend('Full System','Partitioned System');
            set(my_legend, 'Position',[.75 .85 .1 .18]);

        end

        x0 = combine(MIP_cell{i}{1,2},MIP_cell{i}{2,2});
        subtit_pos = pos_BR' - [.08 0 .35 0];
        
        uicontrol('Style', 'text',...
        'String', {mod_mat2str(x0),['phi = ', num2str(min(phi_f,phi_b))]},... %replace something with the text you want
        'Units','normalized',...
        'FontSize',13,...
        'BackgroundColor','w',...
        'Position', subtit_pos);
    
%         subtit_pos = subtit_pos + [.920 0 0 0];
%         
%         uicontrol('Style', 'text',...
%         'String', ['phi = ', num2str(min(phi_f,phi_b))],... %replace something with the text you want
%         'Units','normalized',...
%         'FontSize',10,...
%         'Position', subtit_pos);
    
        set(gcf,'OuterPosition',[200 100 1000 700]);
        
%         figure(fig_st+fig_pi+10)
%         set(gcf,'Position',pos_fig)
%         plot_BRFR(BR_p,FR_p,pos_BR',pos_FR',s_title_p,labelON)
%         if i == 1
%             xlabel(sy);
%             
%         elseif i == N
%             
%             uicontrol('Style', 'text',...
%             'String', 'MIP of Core Concepts For Complex',... %replace something with the text you want
%             'Units','normalized',...
%             'FontSize',16,...
%             'Position', [0.33 0.9 0.43 0.08]);
%         
%             uicontrol('Style', 'text',...
%             'String', 'PAST <-- CURRENT',... %replace something with the text you want
%             'Units','normalized',...
%             'FontSize',13,...
%             'Position', [0.15 0.82 0.28 0.05]);
%         
%             uicontrol('Style', 'text',...
%             'String', 'CURRENT --> FUTURE',... %replace something with the text you want
%             'Units','normalized',...
%             'FontSize',13,...
%             'Position', [0.59 0.82 0.28 0.05]);
%         end
%         
%         x0 = combine(MIP_cell{i}{1,2},MIP_cell{i}{2,2});
%         subtit_pos = pos_BR' - [.08 0 .33 0];
%         
%         uicontrol('Style', 'text',...
%         'String', mod_mat2str(x0),... %replace something with the text you want
%         'Units','normalized',...
%         'FontSize',13,...
%         'Position', subtit_pos);
%     
%         subtit_pos = subtit_pos + [.920 0 0 0];
%         
%         uicontrol('Style', 'text',...
%         'String', ['phi = ', num2str(min(phi_f,phi_b))],... %replace something with the text you want
%         'Units','normalized',...
%         'FontSize',10,...
%         'Position', subtit_pos);
        
        s = [string{3}, ' || ', string_p{3},' || phi = ',num2str(phi_vec(i))];
        fprintf('%s\n',s);
    end
    
    
    
else
    f_size = zeros(1,2);
    if N >= 18
        f_size(1) = 3;
        f_size(2) = 6;
    elseif N >= 9
        f_size(1) = 3;
        f_size(2) = ceil(N/3);
    elseif N < 9
        f_size(1) = 2;
        f_size(2) = ceil(N/2);
    end
    
    r = f_size(1);
    c = f_size(2);
    fig_max = r*c;
    
    save test_set
    
    pos_vec = zeros(fig_max,4);
    m_s = 40;
    p_s = 300;
    margin = [0, m_s, 0, 0 ];
    
    pos_fig = [0,0,p_s*c,p_s*r];
    pos_fig = pos_fig + r*margin;
    
    r_m = m_s/(r*p_s+r*m_s);
    r_p = r*p_s/(r*p_s+r*m_s);
    
    for i=1: r
        for j=1: c
            k = (i-1)*c + j;
            pos_vec(k,:) = [1/c*(j-1), 1/r*(r-i)*r_p, 0, 0];
            pos_vec(k,:) = pos_vec(k,:) + (r-i+1)*[0, r_m, 0, 0];
        end
    end
    
    pos_TR = [0.05/c, 0.15/r*r_p, 0.78/c, 0.78/r*r_p];
    
    pos_BR = r_p*[0.08/c, 0/r*r_p, 0.07/c, 0.005/r*r_p];
    pos_BR = pos_BR + [pos_TR(3),pos_TR(2),0,pos_TR(4)];
    
    pos_FR = r_p*[0.05/c, 0.05/r*r_p, 0.005/c, 0.07/r*r_p];
    pos_FR = pos_FR + [0, 0, pos_TR(3), 0];
    
    pos3 = zeros(3,4);
    
    for i = 1: N
        fig_pi = floor((i-1)/fig_max);
        
        prob = REP_cell{i,1};
        prob_prod = REP_cell{i,2};
        [string_p string] = make_title_fb(MIP_cell{i});
        s_title= [string,': \phi=',num2str(phi_vec(i))];
        
        j = i - fig_pi*fig_max;
        pos_add = pos_vec(j,:);
        
        pos3(1,:) = pos_TR + pos_add;
        pos3(2,:) = pos_BR + pos_add;
        pos3(3,:) = pos_FR + pos_add;
        
        figure(fig_st+fig_pi);
        % figure('Resize','off');
        set(gcf,'Position',pos_fig)
        % set(gcf,'MenuBar','None')
        plot_TRBRFR(prob,pos3,s_title);
        if i == N
            Big_phi_comp = sum(phi_vec);
            sy = ['\Phi=',num2str(Big_phi_comp,3)];
            xlabel(sy)
        end
        
        figure(10+fig_st+fig_pi);
        % figure('Resize','off')
        set(gcf,'Position',pos_fig)
        % set(gcf,'MenuBar','None')
        plot_TRBRFR(prob_prod,pos3,string_p);
        
        if i== 1
            % fprintf('Irreducible points\n');
        end
        fprintf('%s\n',s_title);
       
    end
    
end

end

function [] = plot_BRFR(BR,FR,pos_BR,pos_FR,s_title,labelON)

subplot('Position',pos_BR)
h = bar(0:length(BR)-1,BR,'hist');
set(gca,'XTick',0:length(BR)-1)
set(gca,'YTickLabel',[0 .5 1])

states = convert(length(BR));
axis([-0.5 length(BR)-0.5 0 1.0])
title(s_title{1})
if labelON== 1
    set(gca,'XTickLabel',num2str(states,'%d')) % uncomment this to have a binary valued x-axis
    rotateXLabels( gca(), 90) % uncomment if binary values are used on the x-axis
    xlab = xlabel('State (Node Order: [1...N])','FontSize',14,'Units','pixels');
%     disp(get(xlab,'Position'))
    set(xlab,'Position',get(xlab,'Position') - [0 100 0])
% disp(get(xlab,'Position'))
else
    set(gca,'XTickLabel',[]) 
end


subplot('Position',pos_FR);
h = bar(0:length(FR)-1,FR,'hist');
set(gca,'XTick',0:length(FR)-1)
set(gca,'YTickLabel',[0 .5 1])
states = convert(length(FR));
axis([-0.5 length(FR)-0.5 0 1.0])
title(s_title{2})
if labelON== 1
    set(gca,'XTickLabel',num2str(states,'%d')) % uncomment this to have a binary valued x-axis
    rotateXLabels( gca(), 90) % uncomment if binary values are used on the x-axis
%     xlabel('State (Node Order: [1...N])','FontSize',14)
    xlab = xlabel('State (Node Order: [1...N])','FontSize',14,'Units','pixels');
    set(xlab,'Position',get(xlab,'Position') - [0 100 0])

else
    set(gca,'XTickLabel',[]) 
end


end

function [] = plot_TRBRFR(TR,pos3,s_title)

BR = sum(TR,2);
BR = BR/sum(BR);

FR = sum(TR,1);
FR = FR/sum(FR);

pos_TR = pos3(1,:);
pos_BR = pos3(2,:);
pos_FR = pos3(3,:);

subplot('Position',pos_TR);
imagesc(TR);
colormap(flipud(gray))
caxis([0 1])
set(gca,'XTick',[])
set(gca,'YTick',[])
title(s_title)

subplot('Position',pos_BR)
h = bar(0:length(BR)-1,BR);
set(h,'facecolor','black')
set(gca,'XTick',0:length(BR)-1)
set(gca,'YTickLabel',[])
states = convert(length(BR));
set(gca,'XTickLabel',num2str(states,'%d')) % uncomment this to have a binary valued x-axis
axis([-0.5 length(BR)-0.5 0 1.0])
title('BR')
set(gca,'YAxisLocation','right')
set(gca,'XAxisLocation','top')
set(gca,'CameraUpVector',[-1 0 -1])

subplot('Position',pos_FR);
h = bar(0:length(FR)-1,FR);
set(h,'facecolor','black');
set(gca,'XTick',0:length(FR)-1)
set(gca,'YTickLabel',[])
states = convert(length(FR));
set(gca,'XTickLabel',num2str(states,'%d')) % uncomment this to have a binary valued x-axis
axis([-0.5 length(FR)-0.5 0 1.0])
rotateXLabels( gca(), 90) % uncomment if binary values are used on the x-axis
% title('FR')

end

function states = convert(N)
states = zeros(N,log2(N));
for i=1: N
    sigma = trans2(i-1,log2(N));
    states(i,:) = sigma';
end
end
