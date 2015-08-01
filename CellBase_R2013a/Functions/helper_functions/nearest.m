function pos = nearest(X,y,varargin)
%NEAREST   Find nearest point in vector.
%   POS = NEAREST(X,Y,SIDE) returns the position of an element of vector X
%   that is the nearest to the specified y.
%   If SIDE is specified to be other than zero: 
%   SIDE = -1, from left 
%   SIDE = 1, from right
%
%   See also PLOT_TIMECOURSE and ERRORSHADE.

%   Edit log: BH 6/27/11

% Find nearest
if ~isempty(X)
    if nargin > 2
        tmp = X - y;
        if (varargin{1} > 0)
            tmp = tmp(tmp > 0);
            if ~isempty(tmp)
                pos = find(X-y == min(tmp));
            else
                pos = [];
            end
        elseif varargin{1} < 0   %'left'
            
            tmp = tmp(tmp<=0);
            if ~isempty(tmp)
                pos = find(X-y == max(tmp));
            else
                pos = [];
            end
        else
            pos = find(abs(X-y)==min(abs(X-y)));    % real nearest
        end
    else
        pos = find(abs(X-y)==min(abs(X-y)));    % real nearest
    end
end % isempty

if exist('pos','var')
    pos = min(pos);
else
    pos = [];
end