function cb_name = checknmcb(cb_name)
%CHECKNMCB   Check whether a CellBase name is valid.
%   CB_NAME = CHECKNMCB(CB_NAME) is a helper function. It checks whether
%   CB_NAME is a valid name of an initialized CellBase. It forces to choose
%   a new valid name that is returned, or it exits with an error.
%
%   See also INITCB and DELETECB.

%   Edit log: BH 5/30/11

% Input argument check
error(nargchk(1,1,nargin))

% If this is the first CellBase
if ~ispref('cellbase','cellbases')
    return
end

% Check if the name is used
cellbases = getpref('cellbase','cellbases');
nmcb = length(cellbases);
names = cell(1,nmcb);
for k = 1:nmcb
    names{k} = cellbases{k}(1).name;
end
crcb = find(strcmp(names,cb_name),1);
if ~isempty(crcb)
    disp('CellBase with the specified name already exists.')
    disp('Names of initialized CellBases: ')
    disp(names)
    cb_name = input('Give a unique CellBase name! ','s');
end

% Last check
crcb = find(strcmp(names,cb_name),1);
if ~isempty(crcb)
    error('Name is used.')
end