function  Result = runanalysis(funhandle,varargin)
%RUNANALYSIS   Run an analysis on CellBase.
%   R = RUNANALYSIS(FH,ARGS,'cellids',C) executes an analysis function 
%   (passed by FH, function handle) and returns the result in a matrix (R)
%   with each row being a cell and each column a result of the analysis.
%   Input arguments to FH can be passed in ARGS. The optional pair of input
%   arguments 'cellids', C determines which cell IDs the function will be
%   executed for; by default all cells of the CellBase are included (see
%   CellBase documentation). The 'cellids', C pair can appear at any
%   position of the input argument list.
%
%   R = RUNANALYSIS(FH,ARGS,'cellids',C,'outputargs',OP) accepts the
%   optional input argument pair 'outputargs', OP that determines which
%   output arguments of FH should be included in R (by default, all outputs
%   of FH are included). This input argument pair can appear at any
%   position of the input argument list.
%
%   Examples:
%   runanalysis(@Lratio,{'Energy','Peak'},[1:4],'cellids',CELLIDLIST(27))
%   runanalysis(@Lratio,'cellids',CELLIDLIST(27),{'Energy','Peak'},[1:4],...
%       'outputargs',[1 2])
%
%   See also ADDANALYSIS.

%   Edit log: BH 4/23/2012

% Input argument check
if nargin < 1
	help runanalysis
	return
end
try   % check analysis function
    funname = func2str(funhandle);
catch
    disp('RUNANALYSIS: Function handle is not valid.');
    return
end
if any(strcmp(varargin,'cellids'))   % parse cellids
    cinx = find(strcmp(varargin,'cellids'));
    CELLIDS = varargin{cinx+1};
    arglist = varargin(setdiff(1:length(varargin),[cinx cinx+1]));
else
    CELLIDS = listtag('cells');
    arglist = varargin;
end
if any(strcmp(arglist,'outputargs'))   % parse outputargs
    oinx = find(strcmp(arglist,'outputargs'));
    opargs = arglist{oinx+1};
    arglist = arglist(setdiff(1:length(arglist),[oinx oinx+1]));
else
    opargs = 1:nargout(funhandle);
end

% Load CellBase
load(getpref('cellbase','fname'));

% Run analysis
disp('RUNANALYSIS ...');
NumCells = length(CELLIDS);
Result = cell(NumCells,length(opargs));
for cellnum = 1:NumCells    % cellid loop
    print_progress(cellnum,round(NumCells/100),5);
    
    % Execute analysis
    fullarglist = {char(CELLIDS{cellnum}) arglist{:}};  %#ok<CCAT>
    [property_values{1:nargout(funhandle)}] = feval(funhandle,fullarglist{:});  
    Result(cellnum,:) = property_values(opargs);
end

% Convert to double if possible
if all(cellfun(@(s)isequal(numel(s),1),Result)&cellfun(@isnumeric,Result))
    Result = cell2mat(Result);
end