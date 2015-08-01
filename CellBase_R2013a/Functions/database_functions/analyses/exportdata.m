function  varargout = exportdata(varargin)
%EXPORTDATA   Export all data introduced to CellBase by the 'insertdata' function.
%   EXPORTDATA saves all variables from 'TheMatrix' which corresponds to
%   'insertdata' in ANALYSES to a cell table under the name
%   'CellBase_InsertData' in the CellBase directory (see CellBase
%   documentation).
%
%   EXPORTDATA(FILENAME) exports the data to a file named FILENAME in the
%   current directory.
%
%   DATA = EXPORTDATA('VAR') exports the data to a variable instead of
%   saving.
%
%   See also INSERTDATA.

%   Edit log: BH, 7/6/12

% Load CellBase
load(getpref('cellbase','fname'));

% Find entries corresponding to 'insertdata'
[jnk colpos] = findanalysis(@insertdata);   %#ok<*ASGLU>
[jnk anapos] = findanalysis(@insertdata,'pos');
if ~anapos
    disp('EXPORTDATA: Insertdata has not been used.');
    varargout{1} = [];
    return
end

% Construct output table
NumCell = length(CELLIDLIST);
NumCol = length(colpos);
data = mat2cell(TheMatrix(:,colpos),ones(NumCell,1),ones(NumCol,1)); %#ok<NODEF,MMTC>
propnames = cat(2,ANALYSES(anapos).propnames);
vargin = {ANALYSES(anapos).varargin};
timestamp = {ANALYSES(anapos).timestamp};
vardata{1,1} = 'propnames';
vardata(1,2:length(propnames)+1) = propnames;
vardata{2,1} = 'funhandle';
vardata(2,2:length(propnames)+1) = {@insertdata};
vardata{3,1} = 'varargin';
vardata(3,2:length(vargin)+1) = vargin;
vardata{4,1} = 'timestamp';
vardata(4,2:length(timestamp)+1) = timestamp;
vardata(5:NumCell+4,1) = CELLIDLIST;
vardata(5:NumCell+4,2:NumCol+1) = data;

% Default filename
fname = fullfile(getpref('cellbase','datapath'),'CellBase_InsertData.mat');

% Handle the options specified by the input argument
if nargin == 1
    arg = varargin{1};
    if strncmp(arg,'variable',3)   % return with a variable, no saving
        varargout{1} = vardata;
        return
    elseif ischar(arg)   % redefine filename
        fname = fullfile(pwd,arg);
    else
        disp('EXPORTDATA: Option not understood.');
        return
    end
end

% Save to file
save(fname,'vardata')
disp(sprintf('EXPORTDATA: %d hand entered properties saved in %s',NumCol,fname));