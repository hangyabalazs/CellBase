function pos = restrict2(X,mn,mx,varargin)
%RESTRICT2   Restrict vector to a range. 
%   POS = RESTRICT2(X,MN,MX) returns the positions in vector X that are
%   greater than MN and smaller than MX.
%
%   See also PLOT_TIMECOURSE.

%   Edit log: BH 6/27/11

% Restrict 'X'
if nargin ==  3
    pos0 = find(X>=mn);
    pos1 = find(X(pos0)<=mx);
elseif nargin > 3 
    eval(['pos0=find(X ' varargin{1} ' mn);']);
    if nargin > 4    
        eval(['pos1=find(X(pos0) ' varargin{2} ' mx);']);
    else
        pos1 = find(X(pos0)<=mx);
    end
end
pos = pos0(pos1); 