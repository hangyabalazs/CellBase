function  setmyfigure(figh,varargin)
%SETMYFIGURE   Set figure properties.
%   SETMYFIGURE(F) sets figure properties for F to preset values.
%   Additional parameter settings can be passed as parameter, value pairs
%   (see SET and figure properties).
%
%   See also SETMYPLOT.

%   Edit log: BH 6/27/11

% Set figure properties
set(figh,'Color','w')
set(figh,'DefaultAxesTickDir','out');
set(figh,'DefaultAxesFontName','Ariel','DefaultAxesFontSize',12);
set(figh,'DefaultAxesBox','off')
set(figh,'DefaultAxesColor','none')
set(figh,'DefaultAxesActivePositionProperty','Position');
if nargin > 1
    set(figh,varargin{:},'DefaultAxesFontName','Ariel','DefaultAxesFontSize',12);
end