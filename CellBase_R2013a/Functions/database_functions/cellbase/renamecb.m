function renamecb(cb_name1,cb_name2)
%RENAMECB   Rename CellBase.
%   RENAMECB(NAME1,NAME2) renames CellBase NAME1 to NAME2.
%
%   See also INITCB, CHOOSECB and DELETECB.

%   Edit log: BH 5/31/11

% Check input arguments
error(nargchk(2,2,nargin))

% Find CellBase
crcb = findcb(cb_name1);
cellbases = getpref('cellbase','cellbases');
cb_name1 = cellbases{crcb}.name;

% Rename CellBase
cellbases{crcb}.name = cb_name2;
setpref('cellbase','cellbases',cellbases)

% If the current was renamed
if isequal(getpref('cellbase','name'),cb_name1)
    setpref('cellbase','name',cb_name2)
end