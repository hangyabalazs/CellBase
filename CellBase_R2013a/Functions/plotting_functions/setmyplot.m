function  setmyplot(ax,varargin)
%SETMYPLOT   Set axis properties.
%   SETMYPLOT(AX) sets the font of AX to Ariel (font size, 12) and
%   'TickDir' to 'out'.
%
%   See also SETMYFIGURE.

%   Edit log: BH 6/27/11

% Set axis properties
set(ax,'TickDir','out');
set(ax,'FontName','Ariel','FontSize',12);
if nargin > 1
    set(varargin{:},'FontName','Ariel','FontSize',12);
end