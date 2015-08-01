function cb_name = whichcb
%WHICHCB   Name of the active CellBase.
%   CB_NAME = WHICHCB returns the name of the current CellBase.
%
%   See also CHOOSECB.

%   Edit log: BH 5/30/11

% Get CellBase name
cb_name = getpref('cellbase','name');