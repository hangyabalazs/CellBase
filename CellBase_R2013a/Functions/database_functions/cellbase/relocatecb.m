function relocatecb(cb_name,cb_path)
%RELOCATECB   Change CellBase path.
%   RELOCATECB(NAME,PATH) relocates CellBase NAME to PATH.
%
%   See also INITCB, CHOOSECB, RENAMECB and DELETECB.

%   Edit log: BH 12/31/11

% Check input arguments
error(nargchk(2,2,nargin))

% Find CellBase
crcb = findcb(cb_name);
cellbases = getpref('cellbase','cellbases');
cb_name = cellbases{crcb}.name;
cb_fname = cellbases{crcb}.fname;
[cb_fname_path cb_fname_fname cb_fname_ext] = fileparts(cb_fname);
cb_newfname = fullfile(cb_path,[cb_fname_fname cb_fname_ext]);

% Relocate CellBase
cellbases{crcb}.datapath = cb_path;
cellbases{crcb}.fname = cb_newfname;
setpref('cellbase','cellbases',cellbases)

% If the current was renamed
if isequal(getpref('cellbase','name'),cb_name)
    setpref('cellbase','datapath',cb_path)
    setpref('cellbase','fname',cb_newfname)
end