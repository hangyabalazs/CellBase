function updatecb(cb_name,cb_prop,newvalue)
%UPDATECB   Change CellBase property.
%   UPDATECB(NAME,PROP,VALUE) updates CellBase property PROP to VALUE.
%
%   See also RELOCATECB, INITCB, CHOOSECB, RENAMECB and DELETECB.

%   Edit log: BH 5/22/19

% Check input arguments
narginchk(3,3)

% Find CellBase
crcb = findcb(cb_name);
cellbases = getpref('cellbase','cellbases');
cb_name = cellbases{crcb}.name;

% Update CellBase
cellbases{crcb}.(cb_prop) = newvalue;
setpref('cellbase','cellbases',cellbases)

% If the current was renamed
if isequal(getpref('cellbase','name'),cb_name)
    setpref('cellbase',cb_prop,newvalue)
end