function matches = findcellpos(varargin)
%FINDCELLPOS   Locate animal, session or tetrode positions in CellBase.
%   POS = FINDCELLPOS(CELLID) returns the position of CELLID in
%   CELLIDLIST.
%
%   POS = FINDCELLPOS('RAT',RATID) locates an animal in CELLIDLIST.
%   'ANIMAL' can be used as an alternative for 'RAT' in all sytaxes.
%
%   POS = FINDCELLPOS('RAT',RATID,'SESSION',SESSIONID) locates a session in
%   CELLIDLIST.
%
%   POS = FINDCELLPOS('RAT',RATID,'SESSION',SESSIONID,'TETRODE',TET) 
%   locates a tetrode in CELLIDLIST.
%
%   FINDCELLPOS returns the matching position indices in CellBase for a
%   particular animal, session, cell or tetrode.
%
%   See also FINDALLCELLS and ADDCELL.

%   Edit log: BH 3/21/11, 5/3/12

% Load CellBase
global CELLIDLIST ANALYSES TheMatrix
if isempty(CELLIDLIST)
    load(getpref('cellbase','fname'));
end

% Locate cell ID
if nargin == 1   % most frequent use, requires special handling
    matches = strmatch(varargin{1},CELLIDLIST,'exact');
    if isempty(matches) 
        matches = 0;
    end
    return
end

% Input argument check
N = nargin;
if ~(N==2 || N==4 || N==6)
    argerror;
    return
end

rat = '';
session = '';
tet = -1;
for i = 1:2:N
    switch char(varargin{i})
        case {'rat','animal'}
            rat = char(varargin{i+1});
        case 'session'
            session = char(varargin{i+1});
        case 'tetrode'
            x = varargin{i+1};
            if isnumeric(x)
                tet = x;
            else
                tet = str2num(char(x));
            end
        otherwise
            argerror;
            return
    end
end

% Look for matches
matches = [];
for i = 1:length(CELLIDLIST)
    [ratname,sessionid,tetrode,unit] = cellid2tags(CELLIDLIST{i});
    
    % either there is a match or it's empty & we don't care
    r = strcmp(ratname,rat) | isempty(rat);
    s = strcmp(sessionid,session) | isempty(session);
    t = (tet == -1) | (tetrode==tet); % numeric comparison returns []
    if  r*s*t
        matches = [matches i];
    end
end

if isempty(matches) 
    matches = 0;  
end

% -------------------------------------------------------------------------
function argerror

disp('FINDCELLPOS: There was an error in your arguments.')