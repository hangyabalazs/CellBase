function cellids = findcell(varargin)
%FINDCELL   Locate cells in CellBase.
%   FINDCELL returns the cell ID(s) in CellBase for a particular animal,
%   session or tetrode.
%
%   Syntax:
%   CELLIDS = FINDCELL('RAT',RATID)
%   CELLIDS = FINDCELL('RAT',RATID,'SESSION',SESSIONID)
%   CELLIDS = FINDCELL('RAT',RATID,'SESSION',SESSIONID,'TETRODE',TETRODE)
%
%   See also ADDNEWCELLS, ADDCELL, FINDCELLSTR and FINDALLCELLS.

%   Edit log: BH 10/13/11, 7/6/12

% Load CellBase
global CELLIDLIST ANALYSES TheMatrix
if isempty(CELLIDLIST)
    load(getpref('cellbase','fname'));
end

% Input argument check
nARG = nargin;
iARG = 1;
rat = '';
session = '';
tet = -1;  % comparison of empty numeric is not defined
LIST = '';
while iARG <= nARG; 
    switch lower(char(varargin{iARG}))
        case 'rat'
            if (iARG+1) <= nARG;
                rat = lower(char(varargin{iARG+1}));
            else
                argerror;
            end
            iARG = iARG+2;
        case {'session','ses'}
            if (iARG+1) <= nARG;
                session = lower(char(varargin{iARG+1}));
            else
                argerror;
            end      
            iARG = iARG+2;
        case {'tetrode','tet'}
            if (iARG+1) >= nARG;
                x = varargin{iARG+1};
                if isnumeric(x)
                    tet = x;
                else
                    tet = str2num(char(x));
                end
            else
                argerror;
            end 
            iARG = iARG+2;
        case 'rats'
            LIST = 'rats';
            iARG = iARG + 1;
        case 'sessions'
            LIST = 'sessions';
            iARG = iARG + 1;
        otherwise 
            argerror;
            return
    end  % switch
end  % while nargin

% Output for the 'list' version
NumCells = length(CELLIDLIST);   %#ok<USENS>
biglist = cell(1,NumCells);
if ~isempty(LIST)
    for i = 1:NumCells
        [ratname,sessionid] = cellid2tags(CELLIDLIST{i});
        
        switch LIST
            case 'rats'
                biglist{i} = ratname;
            case 'sessions'
                biglist{i} = strcat(ratname,sessionid);
        end
    end  % i
    
    cellids = unique(biglist);
    return
end

% Find match
matches = [];
for i = 1:length(CELLIDLIST)
    [ratname,sessionid,tetrode] = cellid2tags(CELLIDLIST{i});
    r = strcmpi(ratname,rat) | isempty(rat);   % either there is a match or it's empty
    s = strcmpi(sessionid,session) | isempty(session);
    t = tet==-1 | tetrode==tet;
    if  r * s * t   %#ok<BDLOG>
        matches = [matches i];  %#ok<AGROW>
    end
end
cellids = CELLIDLIST(matches);

% -------------------------------------------------------------------------
function argerror

disp('FINDCELL: There was an error with your arguments.')