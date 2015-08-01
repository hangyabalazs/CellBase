function nout = nargvout(funhandle,funargs)
%NARGVOUT   Number of output arguments.
%   N = NARGVOUT(FUNHANDLE,FUNARGS) returns the number of output arguments
%   for the function FUNHANDLE by calling it with FUNARGS (or the first
%   cell in CELLIDLIST if FUNARGS is empty). Note that NARGOUT is not a
%   valid alternative because here the length of output is determined. 
%
%   See also NARGOUT and ADDANALYSIS.

%   Edit log: BH 3/21/11

% Load CellBase
load(getpref('cellbase','fname'));

% Deal with input arguments to the function of interest
cellid = char(CELLIDLIST(1));
if isempty(funargs)
    arglist = {cellid};
else
    if ~iscell(funargs)
        arglist = {cellid funargs};
    else
        arglist = {cellid funargs{:}};
    end
end

% Call the function
value = feval(funhandle,arglist{:}); 
nout = length(value);