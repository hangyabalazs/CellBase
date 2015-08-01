function  status = setvalue(cellid,property,value)
%SETVALUE   Set value associated with a property.
%   STATUS = SETVALUE(CELLID,PROP,VALUE) sets a property (PROP) to a
%   specified value (VALUE) for a given cell (CELLID) in CellBase. If
%   editing TheMatrix is sucessful, STATUS = true is returned.
%
%   See also GETVALUE.

%   Edit log: BH 4/16/2013

% Load CellBase
global CELLIDLIST ANALYSES TheMatrix
if isempty(CELLIDLIST)
    load(getpref('cellbase','fname'));
end

% Get the position of the searched property
[pos pos0] = findanalysis(property);  %#ok<NASGU>
if pos == 0
    disp('SETVALUE: No matching property found.');
    status = false;
    return
end

% Get the position of the cell
cellpos = findcellpos(cellid);

% Edit TheMatrix
if iscell(TheMatrix(cellpos,pos)) && ~iscell(value)
    value = {value};   % convert to cell if necessary
end
TheMatrix(cellpos,pos) = value;

% Return changed variables to workspace & save all
assignin('base','TheMatrix',TheMatrix)
assignin('base','ANALYSES',ANALYSES)
cb = getpref('cellbase','fname');
[pth fnm ext] = fileparts(cb);
dsr = datestr(now);
dsr = regexprep(dsr,':','_');
backup_name = fullfile(pth,[fnm '_' dsr ext]);
copyfile(cb,backup_name)    % make backup before overwriting
save(cb,'TheMatrix','ANALYSES','CELLIDLIST')

% Feedback
status = true;