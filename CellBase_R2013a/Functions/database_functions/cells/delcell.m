function delcell(cellid)
%DELCELL   Delete cells from CellBase.
%   DELCELL(CELLID) deletes a cell (CELLID, character array) or a lsit of
%   cells (CELLID, cell array of strings) from 'CELLIDLIST' and 'TheMatrix'
%   (see CellBase documentation). CellBase is overwritten after saving a
%   backup file.
%
%   See also ADDCELL, INSERTDATA and DELANALYSIS.

%   Edit log: BH 6/7/12

% Load CellBase
load(getpref('cellbase','fname'));

% Find position(s) in CellBase
if ischar(cellid)
    cellpos = findcellpos(cellid);
elseif iscellstr(cellid)
    NumCells = length(cellid);
    cellpos = nan(1,NumCells);
    for i = 1:NumCells
        cellpos(i) = findcellpos(cellid(i));
    end
else
    disp('ADDCELL: Wrong argument type.');
    return
end
if ~cellpos
    disp('ADDCELL: Cell not found in cellbase.');
    return
end

% Delete from list of cellID's
CELLIDLIST(cellpos) = [];

% Delete rows from 'TheMatrix;
TheMatrix(cellpos,:) = [];

% Return changed variables to workspace & save all
assignin('base','TheMatrix',TheMatrix);
assignin('base','CELLIDLIST',CELLIDLIST);
cb = getpref('cellbase','fname');
[pth fnm ext] = fileparts(cb);
dsr = datestr(now);
dsr = regexprep(dsr,':','_');
backup_name = fullfile(pth,[fnm '_' dsr ext]);
copyfile(cb,backup_name)    % make backup before overwriting
save(cb,'TheMatrix','ANALYSES','CELLIDLIST');
clear global CELLIDLIST ANALYSES TheMatrix

% Feedback
if ischar(cellid)
    donestr = sprintf('DELCELL done.\n Cellid: %s deleted. %d cells remain.\n',cellid,length(CELLIDLIST));
else
    donestr = sprintf('DELCELL done.\n Cellids: %s to %s deleted. %d cells remain.\n',char(cellid(1)),char(cellid(end)),length(CELLIDLIST));
end
disp(donestr);