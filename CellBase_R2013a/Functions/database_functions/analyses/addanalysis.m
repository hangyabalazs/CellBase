function  addanalysis(funhandle,varargin)
%ADDANALYSIS   Add an analysis to CellBase.
%   ADDANALYSIS(FUNHANDLE) adds an analysis to CellBase with the specified
%   function handle. ADDANALYSIS also executes the analysis function passed
%   by FUNHANDLE on the cells included in CellBase. Note that the analysis
%   specified by FUNHANDLE should take a cell ID as its first input
%   argument.
%
%   ADDANALYSIS(FUNHANDLE,'PROPERTY_NAMES',PR) passes the list of property
%   names (PR, cell array of strings) to ADDANALYSIS. The property name
%   'default' is used if no property names are specified. The first N
%   outputs of the analysis function are used, where N is the number of
%   inserted properties (length of PR).
%
%   ADDANALYSIS(FUNHANDLE,'PROPERTY_NAMES',PR,'OUTPUT_SUBSET',OS) passes
%   the list of property names (PR, cell array of strings) to ADDANALYSIS.
%   The additional input argument OS determines which outputs of the
%   analysis function are included. OS should be a numerical array indexing
%   into the output argument list of the analysis function (integers
%   between 1 and the number of FUNHANDLE output arguments - see NARGOUT).
%   The length of OS should match the number of inserted properties (i.e.
%   the length of PR).
%
%   ADDANALYSIS(FUNHANDLE,'PROPERTY_NAMES',PR,'OUTPUT_SUBSET',OS,'MANDATORY',ARGM,'ARGLIST',ARGVAL)
%   also passes input arguments to FUNHANDLE. Mandatory arguments should be 
%   passed in ARGM and parameter-value pairs should be specified in ARGVAL
%   (N-by-2 cell array). These are be passed on to the analysis function.
%
%   Examples:
%   addanalysis(@LRatio2,'property_names',{'ID_PC','Lr_PC'},'arglist',{'fea
%   ture_names' {'WavePC1' 'Energy'}})
%
%   addanalysis(@ultimate_psth,...
%   'mandatory',{'trial' @findAlignEvent_negfeedback_gonogo [-0.6 0.6]},...
%   'property_names',{'FA_psth' 'FA_psth_stats'},'output_subset',[1 6],...
%   'arglist',{'dt',0.001;'display',false;'sigma',0.02;'parts','all';'isadaptive',2;...
%   'event_filter','custom';'filterinput','FalseAlarm==1';'maxtrialno',Inf;...
%   'baselinewin',[-0.5 0];'testwin',[0 0.5];'relative_threshold',0.1});
%
%   See also FINDANALYSIS.

%   Edit log: AK 11/06, BH 3/24/11, 12/6/12, 5/16/13

% Input arguments
prs = inputParser;
addRequired(prs,'funhandle',@(s)isa(s,'function_handle'))   % function handle
addParamValue(prs,'property_names','default',@(s)ischar(s)|iscellstr(s))   % output argument names
addParamValue(prs,'mandatory',{},@iscell)   % mandatory input arguments
addParamValue(prs,'arglist',{},@iscell)   % input arguments parameter-value pairs
addParamValue(prs,'output_subset',[],@isnumeric)   % use a subset of output arguments
parse(prs,funhandle,varargin{:})
g = prs.Results;

% Load CellBase preferences
load(getpref('cellbase','fname'));

% Get default property tags
if strcmpi(g.property_names,'default')
    try
        g.property_names = cellstr(feval(funhandle,'default'));  % default property names
    catch ME
        disp(ME.message)  % no default property names
        error('addanalysis:noDefaultProp','ADDANALYSIS: %s does not provide property names.',func2str(funhandle))
    end
end

% Check number of properties 
funname = func2str(funhandle);   % analysis function name
nout = nargout(funname);  % number of output arguments of the analysis function
Lpropnames = length(g.property_names);   % number of property names
if Lpropnames > nout    % too many property names
    error('addanalysis:tooManyPropNames','ADDANALYSIS: Too many property names.')
end

% Is there an identical analysis already?
[jnk prevanal] = findanalysis(funhandle,'position');  %#ok<ASGLU> % look for analysis in CellBase
if prevanal
    for i = 1:length(prevanal)
        if isequal(ANALYSES(floor(prevanal(i))).varargin,g.arglist)  %#ok<NODEF> % analysis already in CellBase
            error('addanalysis:analysisExists','ADDANALYSIS: Analysis already performed in column %d. Delete first.',prevanal(i))
        end
    end
end

% Is there an identical property name?
for i = 1:Lpropnames
    prevanal = findanalysis(g.property_names{i});  % look for property in CellBase
    if prevanal
        error('addanalysis:propertyExists','ADDANALYSIS: Property name already exists at position %d. Delete first or choose a different one.',prevanal)
    end
end

% Find the position for the new analysis
NumAnal  = length(ANALYSES);  % number of analyses in CellBase
NumCells = length(CELLIDLIST);   %#ok<USENS> % number of cells in CellBase
NewAnal = NumAnal + 1;  % rank of new analysis
if NumAnal == 0
    lastcolumn = 0;
else
    lastcolumn = ANALYSES(NumAnal).columns(end);  % last column in TheMatrix
    [NC, NA] = size(TheMatrix);  %#ok<NODEF> % size of TheMatrix
    if (NC ~= NumCells) || (NA ~= lastcolumn)
        error('addanalysis:databaseError','ADDANALYSIS: Internal database inconsistency!')
    end
end

% Add new analyis
columns = lastcolumn+1:lastcolumn+Lpropnames;  % allocate columns for new properties
pargl = g.arglist';
ANALYSES(NewAnal).funhandle  = funhandle;   % function handle
ANALYSES(NewAnal).varargin   = [g.mandatory pargl(:)'];   % fixed input arguments
ANALYSES(NewAnal).propnames  = g.property_names;   % property names
ANALYSES(NewAnal).output_subset = g.output_subset;   % out arguments to include
ANALYSES(NewAnal).columns    = columns;   % columns in TheMatrix
ANALYSES(NewAnal).timestamp  = timestamp;   % execution time stamp

% Execute analysis
for cellnum = 1:NumCells   % loop through all cells
    disp(CELLIDLIST{cellnum})
    print_progress(cellnum,round(NumCells/100),5);   % progress indicator
    arglist = [char(CELLIDLIST{cellnum}) g.mandatory pargl(:)'];  % input arguments for the analysis
    try
        [property_values{1:nargout(funhandle)}] = feval(funhandle,arglist{:});  % run analysis
    catch ME
        disp(ME.message)
        disp('ADDANALYSIS: Error while executing the analysis function. Adding NaNs.')
        property_values = num2cell(nan(1,nout));   % if there is a error in execution, insert a NaN
    end
    
    % Insert into TheMatrix
    if iscell(TheMatrix)
        if ~isempty(g.output_subset)
            TheMatrix(cellnum,columns) = property_values(g.output_subset);
        else
            TheMatrix(cellnum,columns) = property_values(1:Lpropnames);
        end
    else
        if ~isempty(g.output_subset)
            TheMatrix(cellnum,columns) = cell2mat(property_values(g.output_subset));
        else
            TheMatrix(cellnum,columns) = cell2mat(property_values(1:Lpropnames));
        end
    end
end

% Feedback
if Lpropnames == 1
    noutSTR = 'y';
else
    noutSTR = 'ies';
end
donestr = sprintf('\nADDANALYSIS done.\n  %s executed on %d cells creating %d new propert%s.\n',func2str(funhandle),cellnum,Lpropnames,noutSTR);
disp(donestr);

% Return changed variables to workspace & save all
assignin('base','TheMatrix',TheMatrix)
assignin('base','ANALYSES',ANALYSES)
cb = getpref('cellbase','fname');
[pth fnm ext] = fileparts(cb);
dsr = datestr(now);
dsr = regexprep(dsr,':','_');
backup_name = fullfile(pth,[fnm '_' dsr ext]);  % time stamped backup
copyfile(cb,backup_name)    % make backup before overwriting

% SAVE CELLBASE
save(getpref('cellbase','fname'),'TheMatrix','ANALYSES','CELLIDLIST')
clear global CELLIDLIST ANALYSES TheMatrix

% -------------------------------------------------------------------------
function  ts = timestamp

c = clock;
ts = sprintf('%d/%d/%d %.2d:%.2d',c(2),c(3),c(1),c(4),c(5));

% -------------------------------------------------------------------------
function  out = cell2num(arg)

try
    out = cell2mat(arg);
    if ~isnumeric(out)
        out = arg;
    end
catch
    out = arg;
end