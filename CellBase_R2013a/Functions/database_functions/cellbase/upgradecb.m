function upgradecb
%UPGRADECB   Upgrade CellBase to new version.
%   UPGRADECB creates new CellBase preferences to render old instantiations
%   of CellBase compatible with new CellBase functions.
%
%   See also INITCB.

%   Edit log: BH 8/20/13

% Timestamp conversion
tsc = questdlg('I''m going convert timestamps to seconds in','Timestamp conversion',...
    'CellBase','MClust','CellBase');
switch tsc
    case 'CellBase'
        timefactor = 1e-4;  % TT*.mat files will reflect the timestamps of the Ntt files;
                            % loadcb will convert them to seconds
    case 'MClust'
        timefactor = 1;     % timestamps in TT*.mat files will already be converted to seconds
end

% Name CellBase
cb_name = input('Give a name to your CellBase! ','s');
cb_name = checknmcb(cb_name);

% Store cellbases to allow multiple instances
if ispref('cellbase','cellbases')
    cellbases = getpref('cellbase','cellbases');
else
    cellbases = {};
end
gp = getpref('cellbase');
if isfield(gp,'cellbases')
    gp = rmfield(gp,'cellbases');
end
cellbases{end+1} = gp;

% Set preferences (persistent and maintain their values between MATLAB sessions)
setpref('cellbase','timefactor',timefactor);
setpref('cellbase','cellbases',cellbases);
setpref('cellbase','name',cb_name);