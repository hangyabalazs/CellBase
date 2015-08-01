function  delanalysis(funhandle)
%DELANALYSIS   Delete an analysis from CellBase.
%   DELANALYSIS(FUNHANDLE) deletes a particular analysis (passed by
%   funhandle) from CellBase and the column(s) associated with it. If
%   multiple analyses found, it deletes the first.
%
%   DELANALYSIS(PROPERTYNAME) deletes a property in a similar fashion
%   (PROPERTYNAME should be a character array).
%
%   See also ADDANALYSIS and FINDANALYSIS.

%   Edit log: BH 6/21/12

% Load CellBase
load(getpref('cellbase','fname'));

% Delete property
if ischar(funhandle)  % property
    propname = funhandle;
    [pos pnum] = findanalysis(propname); %#ok<*ASGLU>
    if pnum ~= 0
        if length(ANALYSES(pnum(1)).propnames) > 1 %#ok<NODEF>
            disp(['DELANALYSIS: Multiple properties belong to the same analysis: ' func2str(ANALYSES(pnum).funhandle)]);
            disp('Property cannot be deleted by itself, need to delete entire analysis.');
            return
        end
    end
    anum = pnum;
else  % if not property
    
    % Find analysis
    [pos anum] = findanalysis(funhandle,'position');
end

% Check if analysis was found
if anum == 0
    disp('DELANALYSIS: Analysis function / property not found.');
    return
end
if length(anum) > 1
    anum = anum(1);
    disp('DELANALYSIS: multiple analyses functions found. deleting the first.');
end
cols = ANALYSES(anum).columns;

% Delete from ANALYSES & change column ordering
ANALYSES(anum) = [];
if anum <= length(ANALYSES)    %if there are more analyses columns need to be renumbered
    diffcol = ANALYSES(anum).columns(1) - cols(1);
    for i = anum:length(ANALYSES)
        ANALYSES(i).columns =  ANALYSES(i).columns - diffcol; %#ok<AGROW>
    end
end

% Delete appropriate columns from TheMatrix
TheMatrix(:,cols) = [];

% Output message
if ~ischar(funhandle)
    xstr = func2str(funhandle);
else
    xstr = funhandle;
end
donestr = sprintf('DELANALYSIS done.\n %s deleted from %d cells erasing %d properti(es).\n',xstr,length(CELLIDLIST),length(cols));
disp(donestr);

% Save CellBase
assignin('base','TheMatrix',TheMatrix)
assignin('base','ANALYSES',ANALYSES)
cb = getpref('cellbase','fname');
[pth fnm ext] = fileparts(cb);
dsr = datestr(now);
dsr = regexprep(dsr,':','_');
backup_name = fullfile(pth,[fnm '_' dsr ext]);
copyfile(cb,backup_name)    % make backup before overwriting
save(getpref('cellbase','fname'),'TheMatrix','ANALYSES','CELLIDLIST') 
clear global CELLIDLIST ANALYSES TheMatrix