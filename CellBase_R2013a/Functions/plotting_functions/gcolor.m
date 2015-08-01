function color_str = gcolor(number,varargin)
%GCOLOR   Get color string.
%   S = GCOLOR(N) returns a color string according to a predefined color
%   set.
%
%   S = GCOLOR(N,K) returns a color string and a marker symbol according to
%   predefined color and marker set. Either of the input numbers can be
%   replaced with the corresponding color string or marker symbol.
%
%   See also PLOT.

%   Edit log: BH 1/28/13

% Define color and marker sets
colorV = ['y';'b';'r';'g';'k';'c';'m'];
markerV = ['v','s','o','^','d','h','<','p','>'];

% Input argument handling
if nargin == 1
    type = '';
else
    type = varargin{1};
    if ~ischar(type)
        typeNUM = mod(type,length(markerV));
        type = markerV(typeNUM);   % marker string
    end
end

% Color
if ischar(number)
    [COLOR_FOUND, loc] = ismember(number,colorV);   % find color
    if COLOR_FOUND
        colorX = loc - 1;   % 1 will be added later to avoid 0 from MOD
    else
        disp('Error in GCOLOR: color not found.')
    end
else
    colorX = mod(number,length(colorV));   % color string
end

% Output string
color_str = strcat(colorV(colorX+1),type);   % color + marker