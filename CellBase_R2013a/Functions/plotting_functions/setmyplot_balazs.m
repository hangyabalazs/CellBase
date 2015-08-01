function  setmyplot_balazs(ax,varargin)
%SETMYPLOT_BALAZS   Set axis properties.
%   SETMYPLOT_BALAZS(AX) sets axis properties.
%
%   See also SETMYPLOT and SETMYFIGURE.

%   Edit log: BH 9/29/13

% Input arguments
if nargin == 0
    ax = gca;
end

% Set axis properties
set(ax,'TickDir','out','box','off');
set(ax,'FontSize',12,'LineWidth',0.75);
set(ax,'TickLength',[0.025 0])
set(ax,'FontName','Arial');

% Title and axis label font size
th = findobj(allchild(gca),'Type','text','VerticalAlignment','bottom','Rotation',0);
axlab = setdiff(findobj(allchild(gca),'Type','text'),th);
if ~isempty(th)
    set(th,'FontSize',24,'FontName','Arial')
end
if ~isempty(axlab),
    set(axlab,'FontSize',13,'FontName','Arial')
end
if nargin > 1
    set(varargin{:});
end