function crcb = findcb(cb_name)
%FINDCB   Find CellBase.
%   I = FINDCB(CB_NAME) is a helper function. It finds CellBase with name
%   CB_NAME in the list of initialized instances of CellBase and returns
%   its index. It forces to choose a valid name or returns with an error. 
%
%   See also WHOSCB, CHOOSECB, DELETECB and RENAMECB.

%   Edit log: BH 5/30/11

% Input argument check
error(nargchk(0,1,nargin))
if nargin < 1
    cb_name = '';
end

% Find CellBase
cellbases = getpref('cellbase','cellbases');
nmcb = length(cellbases);
names = cell(1,nmcb);
for k = 1:nmcb
    names{k} = cellbases{k}(1).name;
end
if isempty(cb_name)
    disp('Names of initialized CellBases: ')
    disp(names)
    cb_name = input('Give a valid CellBase name! ','s');
end

% Check for match
crcb = find(strcmp(names,cb_name),1);
if isempty(crcb)
    disp('No match with initialized CellBases.')
    disp('Names of initialized CellBases: ')
    disp(names)
    cb_name = input('Give a valid CellBase name! ','s');
end

% Last check
crcb = find(strcmp(names,cb_name));
if isempty(crcb)
    error('No match with initialized CellBases.')   % stupid case
end