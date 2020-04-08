function OK = addcell(cellid,varargin)
%ADDCELL   Add a cell to the database and perform all the analyses.
%   OK = ADDCELL(CELLID) adds the specified cell to the database and
%   performs the analyses. CELLID may also be a cell array of strings.
%   ADDCELL returns 1 if no errors occur and 0 otherwise.
%
%   OK = ADDCELL(CELLID,'QUIET') suppresses text displays.
%
%   See also ADDNEWCELLS.

%   Edit log: AK 8/05, BH 3/21/11, 5/20/13, 4/18/14

% Input argument check
VERBOSE = 1;
if nargin > 1
    if strcmpi(varargin{1},'quiet')
        VERBOSE = 0;
    end
end

% For multiple cellids
OK = 0;
if iscellstr(cellid)
    for i = 1:length(cellid)
        addcell(char(cellid(i)));
    end
    return
end

% Load cellbase
clear global CELLIDLIST ANALYSES TheMatrix
cellbase_fname = getpref('cellbase','fname');
load(cellbase_fname);

% Determine if new
if findcellpos(cellid)
    disp(sprintf('ADDCELL: Cell %s already in cellbase.',cellid));
    return
end
global NEWCELLIDS
% NEWCELLIDS = [NEWCELLIDS {cellid}];   % enables adding analysis (with addanalysis) that runs iscellid

% Perform analyses
NumCells = length(CELLIDLIST);
NumAnal  = length(ANALYSES);
NewCell  = NumCells + 1;
for i = 1: NumAnal
    
    funhandle = ANALYSES(i).funhandle;
    columns   = ANALYSES(i).columns;
    varg      = ANALYSES(i).varargin;
    
    clear property_values
    if strcmpi(func2str(funhandle),'insertdata')   % we have to handle this separately
        property_name = char(ANALYSES(i).propnames);
%         property_values = {input(sprintf('Enter value for %s ! ',property_name))};
        property_values = num2cell(nan(length(columns),1));
        if size(TheMatrix,1) > 1 && any(cellfun(@ischar,TheMatrix(1,columns)))
            charinx = cellfun(@ischar,TheMatrix(1,columns));   % for character type properties, initialize with empty matrix
            property_values(charinx) = {''};
        end
    else
        if ~isempty(varg)
            try
                if ~isempty(varg)
                    [property_values{1:nargout(funhandle)}] = feval(funhandle,cellid,varg{:});
                else
                    [property_values{1:nargout(funhandle)}] = feval(funhandle,cellid);
                end
            catch
%                 property_values = {nan(length(columns),1)};
                property_values = num2cell(nan(nargout(funhandle),1));
                if size(TheMatrix,1) > 1 && any(cellfun(@ischar,num2cell(TheMatrix(1,columns))))
                    charinx = cellfun(@ischar,TheMatrix(1,columns));   % for character type properties, initialize with empty matrix
                    property_values(charinx) = {''};
                end
            end
        else
            try
                [property_values{1:nargout(funhandle)}] = feval(funhandle,cellid);
            catch
                property_values = num2cell(nan(length(columns),1));
                if iscell(TheMatrix) && size(TheMatrix,1) > 1 && any(cellfun(@ischar,TheMatrix(1,columns)))
                    charinx = cellfun(@ischar,TheMatrix(1,columns));   % for character type properties, initialize with empty matrix
                    property_values(charinx) = {''};
                end
            end
        end
    end
    property_values = property_values(:);  % vector standardized
    if iscell(TheMatrix)   % insert into TheMatrix
        if isfield(ANALYSES,'output_subset') && ~isempty(ANALYSES(i).output_subset) ...
                && ~all(cellfun(@isnan,property_values))
            TheMatrix(NewCell,columns) = property_values(ANALYSES(i).output_subset);
        else
            TheMatrix(NewCell,columns) = property_values(1:length(columns));
        end
    else
        
        if isfield(ANALYSES,'output_subset') && ~isempty(ANALYSES(i).output_subset) ...
                && ~all(cellfun(@isnan,property_values))
            TheMatrix(NewCell,columns) = cell2mat(property_values(ANALYSES(i).output_subset));
        else
            TheMatrix(NewCell,columns) = cell2mat(property_values(1:length(columns)));
        end
    end
end
if NumAnal == 0
    columns = 0;
end
if VERBOSE
    donestr = sprintf('ADDCELL done.\n Cellid ''%s'' added to the database with %d analyses executed, creating %d properti(es).\n',cellid,NumAnal,columns(end));
    disp(donestr);
end
CELLIDLIST{NewCell} = cellid;

% Save cellbase
assignin('base','TheMatrix',TheMatrix)
assignin('base','CELLIDLIST',CELLIDLIST)
save(cellbase_fname,'TheMatrix','ANALYSES','CELLIDLIST')
OK = 1;