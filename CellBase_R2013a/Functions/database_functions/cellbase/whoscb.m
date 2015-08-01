function names = whoscb
%WHOSCB   List of CellBases.
%   L = WHOSCB returns the list of initialized instances of CellBase. 
%
%   See also FINDCB.

%   Edit log: BH 5/31/11

% Input argument check
error(nargchk(0,0,nargin))

% List of CellBase names
cellbases = getpref('cellbase','cellbases');
nmcb = length(cellbases);
names = cell(1,nmcb);
for k = 1:nmcb
    names{k} = cellbases{k}(1).name;
end
disp('Names of initialized CellBases: ')
disp(names)

% Supress output
if nargout < 1
    clear names
end