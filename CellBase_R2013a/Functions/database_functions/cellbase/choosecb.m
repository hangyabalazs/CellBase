function choosecb(cb_name)
%CHOOSECB   Choose between CellBases.
%   CHOOSECB can be used to choose between initialized instances of
%   CellBase.
%
%   CHOOSECB(NAME) sets CellBase to the one with the given NAME.
%
%   See also SWITCHCB and WHICHCB.

%   Edit log: BH 5/30/11

% Input argument check
error(nargchk(0,1,nargin))
if nargin < 1
    cb_name = '';
end

% Find CellBase
crcb = findcb(cb_name);
cellbases = getpref('cellbase','cellbases');

% Set preferences
fld = fieldnames(cellbases{crcb});
fld = setdiff(fld,'cellbases');
for k = 1:length(fld)
    setpref('cellbase',fld{k},cellbases{crcb}(1).(fld{k}))
end

% Set globals
global CELLIDLIST ANALYSES TheMatrix %#ok<NUSED>
if ~isempty(CELLIDLIST)
    fprintf('Updating global variables............')
    loadcb
    fprintf('\b.      Done \n')
end