function [ax, height, extra_plots] = conceptscatter(x,nWholeConcepts, w_highlight_indices, parent_panel, mip_axes)
% BASED ON GPLOTMATRIX

% x = rand(size(in_data))
% assignin('base','x',x);


num_dims = min(size(x,2),8);
num_nodes = log2(num_dims);
concept_var = var(x);
% concept_var_states = zeros(num_dims,1);

replot = ~isempty(mip_axes);

[sorted_vars concept_var_states] = sort(concept_var,'descend');

whole = x(1:nWholeConcepts,:);
part = x(nWholeConcepts+1:end,:);

% whole = rand(size(x(1:nWholeConcepts,:)));
% part = rand(size(x(nWholeConcepts+1:end,:)));
dims = size(x,2);
% 
% assignin('base','whole',whole)
% assignin('base','part',part)

% rows = size(x,2); cols = rows;
rows = 8; cols = rows;
extra_plots = rows - dims;
XvsX = true;


% [g,gn] = mgrp2idx(g,size(x,1),',');
% ng = max(g);

% ynam = xnam;

% Create/find BigAx and make it invisible
% clf; % will this clear the whole gui? YES!
% BigAx = handles.overview_axes;
hold_state = ishold;
% set(BigAx,'Visible','off','color','none','Parent',parent_panel)

clr = 'bgrcmyk';
sym = '.';


siz = repmat(get(0,'defaultlinemarkersize'), size(sym));
% % % if any(sym=='.'),
% % %   units = get(BigAx,'units');
% % % %   set(BigAx,'units','pixels');
% % %   pos = get(BigAx,'Position');
% % % %   set(BigAx,'units',units);
% % %   siz(sym == '.') = max(1,min(15, ...
% % %                    round(15*min(pos(3:4))/size(x,1)/max(rows,cols))));
% % % end


% Store global data for datatips into BixAx
% ginds = cell(1,ng);
% for i=1:ng
%     ginds{i} = find(g==i);
% end

% setappdata(BigAx,'ginds',ginds);
% setappdata(BigAx,'xnam',xnam);
% setappdata(BigAx,'ynam',ynam);
% % % setappdata(BigAx,'x',x);
% % % setappdata(BigAx,'y',x);
% % % setappdata(BigAx,'XvsX',XvsX);
% setappdata(BigAx,'gn',gn);

% TOOK OUT TO GET LINKING WORKING... NOT SURE IF IT MAKES A DIFF
% Make datatips show up in front of axes
% dcm_obj = datacursormode(ancestor(BigAx,'figure'));
% set(dcm_obj,'EnableAxesStacking',true);
% 
% dataCursorBehaviorObj = hgbehaviorfactory('DataCursor');
% set(dataCursorBehaviorObj,'UpdateFcn',@gplotmatrixDatatipCallback);

% Create and plot into axes
ax2filled = false(rows,1);
% % % % 
pos = get(parent_panel,'Position');
width = pos(3)/8;
% height = pos(4)/rows;
height = width;
space = .04; % 2 percent space between axes
pos(1:2) = pos(1:2) + space*[width height]; % shift starting point by spacing
[m,n,k] = size(x); %#ok<ASGLU>
xlim = repmat(cat(3,zeros(rows,1),ones(rows,1)),[rows 1 1]);
ylim = repmat(cat(3,zeros(rows,1)',ones(rows,1)'),[1 cols 1]);


x_bound = [0 1 0];
y_bound = [1 0 0];
% these are the loops that need to be changed to enable data linking

% ax = cell(nchoosek(size(x,2),2)+1,1); % all pairs of dims plus the 3D plot

if replot
    ax = mip_axes;
else
    ax = cell(nchoosek(num_dims,2)+1,1); % all pairs of dims plus the 3D plot
end

ax_index = 1;
for i = 8:-1:0 % count down from rows to 1
   for j = i-1:-1:1, % count down from cols to 1
       
       if ~replot
            axPos = [(j-1)*width+space (rows-i)*height+space ...
                width*(1-space) height*(1-space)];
            ax{ax_index} = axes('Position',axPos, 'visible', 'on', 'Box','on','Parent',parent_panel,...
                'DrawMode','fast','Clipping','On');
            
            xlim(i,j,:) = get(ax{ax_index},'xlim');
            ylim(i,j,:) = get(ax{ax_index},'ylim');
        else
           plots = findobj(ax{ax_index},'Parent',ax{ax_index});
           for k = 1:length(plots)
               delete(plots(k))
           end
        end

        
        if (i <= num_dims)
            
           set(ax{ax_index},'Visible','on')
           state1 = concept_var_states(i);
           state2 = concept_var_states(j);

            plot(ax{ax_index},whole(:,state2),...
                whole(:,state1),'dg')
            



            hold on;
            plot(ax{ax_index},part(:,state2), ...
                part(:,state1),'xb');
            hold on;
            
            plot(ax{ax_index},whole(w_highlight_indices,state2), ...
                whole(w_highlight_indices,state1),'.r');
            hold on;
            
            choices = nchoosek([1 2 3],2);

            for k = 1:size(choices,1)

                hold on
                plot(ax{ax_index},x_bound(choices(k,:)),y_bound(choices(k,:)),'k');

            end
        else
            set(ax{ax_index},'Visible','off')
        end


        
        if j == 1 && i <= num_dims
            ylabel(ax{ax_index},dec2bin(state1-1,num_nodes))
        end
        if i == num_dims && j <= num_dims
            xlabel(ax{ax_index},dec2bin(state2-1,num_nodes))
        end
       
        set(ax{ax_index},'xlimmode','manual','ylimmode','manual','xgrid','off','ygrid','off',...
                'xlim',[-.25 1.25],'ylim',[-.25 1.25],'xticklabel','','yticklabel','')
        ax_index = ax_index + 1;

   end
end

j = ceil(rows/2);

if ~replot
    axPos = [(j+1)*width+space (rows - j - .5)*height+space ...
            width*(1-space)*(j+.75) height*(1-space)*(j+.75)];
    axes3D = axes('Position',axPos, 'visible', 'on', 'Box','on','Parent',parent_panel,'DrawMode','fast');

    ax{ax_index} = axes3D;
else
   plots = findobj(ax{ax_index},'Parent',ax{ax_index});
   for k = 1:length(plots)
       delete(plots(k))
   end
end
scatter3(ax{ax_index},whole(:,concept_var_states(1)),whole(:,concept_var_states(2)),...
    whole(:,concept_var_states(3)),'Marker','d','MarkerFaceColor','g')

hold on

scatter3(ax{ax_index},part(:,concept_var_states(1)),part(:,concept_var_states(2)),...
    part(:,concept_var_states(3)),'Marker','x','MarkerFaceColor','b')

hold on

scatter3(ax{ax_index},whole(w_highlight_indices,concept_var_states(1)),whole(w_highlight_indices,concept_var_states(2)),...
    whole(w_highlight_indices,concept_var_states(3)),'Marker','.','MarkerFaceColor','r')

xlabel(ax{ax_index},dec2bin(concept_var_states(1)-1,num_nodes))
ylabel(ax{ax_index},dec2bin(concept_var_states(2)-1,num_nodes))
zlabel(ax{ax_index},dec2bin(concept_var_states(3)-1,num_nodes))


set(ax{ax_index},'xlimmode','manual','ylimmode','manual',...
        'xlim',[-.25 1.25],'ylim',[-.25 1.25],'zlim',[-.25 1.25],...
        'CameraViewAngleMode','manual')
        


% plot tetrahedron bounds
x_bound = [0 0 1 0];
y_bound = [0 1 0 0];
z_bound = [0 0 0 1];
choices = nchoosek([1 2 3 4],2);

for i = 1:size(choices,1)
    
    hold on
    plot3(ax{ax_index},x_bound(choices(i,:)),y_bound(choices(i,:)),z_bound(choices(i,:)),'k');
    
end

% linkdata on

% replace with real data
% whole = x(1:nWholeConcepts,:);
% part = x(nWholeConcepts+1:end,:);
% assignin('base','whole',whole)
% assignin('base','part',part)


% x(:) = in_data(:);
% assignin('base','x',x);


% ld = linkdata(gcf)
% fieldnames(ld)

% xlimmin = min(xlim(:,:,1),[],1); xlimmax = max(xlim(:,:,2),[],1);
% ylimmin = min(ylim(:,:,1),[],2); ylimmax = max(ylim(:,:,2),[],2);
% 
% % % Set all the limits of a row or column to be the same and leave 
% % % just a 5% gap between data and axes.
% inset = .05;
% for i=2:rows,
%   set(ax(i,1),'ylim',[ylimmin(i,1) ylimmax(i,1)])
%   dy = diff(get(ax(i,1),'ylim'))*inset;
%   set(ax(i,1:i-1),'ylim',[ylimmin(i,1)-dy ylimmax(i,1)+dy])
% end
% for j=1:cols-1,
%   set(ax(j+1,j),'xlim',[xlimmin(1,j) xlimmax(1,j)])
%   dx = diff(get(ax(1,j),'xlim'))*inset;
%   set(ax(j+1:rows,j),'xlim',[xlimmin(1,j)-dx xlimmax(1,j)+dx])
%   if ax2filled(j)
%      set(ax2(j),'xlim',[xlimmin(1,j)-dx xlimmax(1,j)+dx])
%   end
% end



% % Label plots one way or the other
% if (donames && ~isempty(xnam))
%    for j=1:cols
%       set(gcf,'CurrentAx',ax(j,j));
%       h = text((...
%           xlimmin(1,j)+xlimmax(1,j))/2, (ylimmin(j,1)+ylimmax(j,1))/2, -.1,...
%           xnam{j}, 'HorizontalAlignment','center',...
%           'VerticalAlignment','middle');
%    end
% else
%    if ~isempty(xnam)
%       for j=1:cols, xlabel(ax(rows,j),xnam{j}); end
%    end
%    if ~isempty(ynam)
%       for i=1:rows, ylabel(ax(i,1),ynam{i}); end
%    end
% end

% Ticks and labels on outer plots only
% set(ax(1:rows-1,:),'xticklabel','')
% set(ax(:,2:cols),'yticklabel','')
% set(BigAx,'XTick',get(ax(rows,1),'xtick'),'YTick',get(ax(rows,1),'ytick'), ...
%           'userdata',ax,'tag','PlotMatrixBigAx')

% Create legend if requested; base it on the top right plot
% if (doleg)
%    gn = gn(ismember(1:size(gn,1),g),:);
%    legend(ax(1,cols),gn);
% end

% Make BigAx the CurrentAxes
% set(gcf,'CurrentAx',BigAx)
% if ~hold_state,
%    set(gcf,'NextPlot','replace')
% end




% Also set Title and X/YLabel visibility to on and strings to empty
% set([get(BigAx,'Title'); get(BigAx,'XLabel'); get(BigAx,'YLabel')], ...
%  'String','','Visible','on')

% if nargout~=0,
%   h = hh;
%   if any(ax2filled)
%      ax = [ax; ax2(:)'];
%   end
% end

% -----------------------------
function datatipTxt = gplotmatrixDatatipCallback(obj,evt)

target = get(evt,'Target');
ind = get(evt,'DataIndex');
pos = get(evt,'Position');

dtcallbackdata = getappdata(target,'dtcallbackdata');
[BigAx,gnum,row,col] = dtcallbackdata{:};

ginds = getappdata(BigAx,'ginds');
xnam = getappdata(BigAx,'xnam');
ynam = getappdata(BigAx,'ynam');
xdat = getappdata(BigAx,'x');
ydat = getappdata(BigAx,'y');
XvsX = getappdata(BigAx,'XvsX');
gn = getappdata(BigAx,'gn');

gind = ginds{gnum};
obsind = gind(ind);

xvals = xdat(obsind,:);
yvals = ydat(obsind,:);

x = xvals(col);
y = yvals(row);

if x~=pos(1) || y~=pos(2)
    % Something is inconsistent, display default datatip.
    datatipTxt = {sprintf('X: %s',num2str(pos(1))),sprintf('Y: %s',num2str(pos(2)))};
else
    if isempty(xnam)
        xnam = cell(size(xdat,2),1);
        for i = 1:size(xdat,2)
            xnam{i} = sprintf('xvar%s',num2str(i));
        end
    end
    if isempty(ynam)
        ynam = cell(size(ydat,2),1);
        for i = 1:size(ydat,2)
            ynam{i} = sprintf('yvar%s',num2str(i));
        end
    end

    % Generate datatip text.
    datatipTxt = {
        [xnam{col},': ',num2str(x)],...
        [ynam{row},': ',num2str(y)],...
        '',...
        sprintf('Observation: %s',num2str(obsind)),...
        };

    if ~isempty(gn)
        datatipTxt{end+1} = ['Group: ',gn{gnum}];
    end
    datatipTxt{end+1} = '';

    xnamTxt = cell(length(xvals),1);
    for i=1:length(xvals)
        xnamTxt{i} = [xnam{i} ': ' num2str(xvals(i))];
    end
    datatipTxt = {datatipTxt{:}, xnamTxt{:}};
    
    if ~XvsX
        ynamTxt = cell(length(yvals),1);
        for i=1:length(yvals)
            ynamTxt{i} = [ynam{i} ': ' num2str(yvals(i))];
        end
        datatipTxt = {datatipTxt{:}, ynamTxt{:}};
    end

end

function [ogroup,glabel,gname,multigroup] = mgrp2idx(group,rows,sep); 
%MGRP2IDX Convert multiple grouping variables to index vector 
%   [OGROUP,GLABEL,GNAME,MULTIGROUP] = MGRP2IDX(GROUP,ROWS) takes 
%   the inputs GROUP, ROWS, and SEP.  GROUP is a grouping variable (numeric 
%   vector, string matrix, or cell array of strings) or a cell array 
%   of grouping variables.  ROWS is the number of observations. 
%   SEP is a separator for the grouping variable values. 
% 
%   The output OGROUP is a vector of group indices.  GLABEL is a cell 
%   array of group labels, each label consisting of the values of the 
%   various grouping variables separated by the characters in SEP. 
%   GNAME is a cell array containing one column per grouping variable 
%   and one row for each distinct combination of grouping variable 
%   values.  MULTIGROUP is 1 if there are multiple grouping variables 
%   or 0 if there are not. 
 
%   Tom Lane, 12-17-99 
%   Copyright 1993-2002 The MathWorks, Inc.  
%   $Revision: 1.4 $  $Date: 2002/02/04 19:25:44 $ 
 
multigroup = (iscell(group) & size(group,1)==1); 
if (~multigroup) 
   [ogroup,gname] = grp2idx(group); 
   glabel = gname; 
else 
   % Group according to each distinct combination of grouping variables 
   ngrps = size(group,2); 
   grpmat = zeros(rows,ngrps); 
   namemat = cell(1,ngrps); 
    
   % Get integer codes and names for each grouping variable 
   for j=1:ngrps 
      [g,gn] = grp2idx(group{1,j}); 
      grpmat(:,j) = g; 
      namemat{1,j} = gn; 
   end 
    
   % Find all unique combinations 
   [urows,ui,uj] = unique(grpmat,'rows'); 
    
   % Create a cell array, one col for each grouping variable value 
   % and one row for each observation 
   ogroup = uj; 
   gname = cell(size(urows)); 
   for j=1:ngrps 
      gn = namemat{1,j}; 
      gname(:,j) = gn(urows(:,j)); 
   end 
    
   % Create another cell array of multi-line texts to use as labels 
   glabel = cell(size(gname,1),1); 
   if (nargin > 2) 
      nl = sprintf(sep); 
   else 
      nl = sprintf('\n'); 
   end 
   fmt = sprintf('%%s%s',nl); 
   lnl = length(fmt)-3;        % one less than the length of nl 
   for j=1:length(glabel) 
      gn = sprintf(fmt, gname{j,:}); 
      gn(end-lnl:end) = []; 
      glabel{j,1} = gn; 
   end 
end 

