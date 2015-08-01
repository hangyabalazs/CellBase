function  status = writedata(fname)
%WRITEDATA   Write CellBase to Excel.
%   STATUS = WRITEDATA writes CELLIDLIST and TheMatrix (see CellBase
%   documentation) to Excel. STATUS is true if writing Excel file was
%   successful and false otherwise.
%
%   STATUS = WRITEDATA(FNAME) writes the file with the specified file name
%   (FNAME).
%
%   See also INSERTDATA and EXPORTDATA.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   12-Apr-2013

%   Edit log: BH, 4/12/13

% Output file name
if nargin < 1
    [filename, pathname] = uiputfile( ...
        {'*.xls;*.xlsx','Supported formats (*.xls,*.xlsx)';...
        '*.xls;*.xlsx','Excel files (*.xls,*.xlsx)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Export as');     % select folder and file name
    if isequal(filename,0)   % action cancelled
        disp('Data export cancelled.')
        status = false;  % no outfile
        return
    else
        fname = fullfile(pathname,filename);   % output file name   
    end
end
[pt ft ext] = fileparts(fname);
if isempty(ext)
    ext = '.xls';   % default output file format
    fname = fullfile(pt,[ft ext]);   % output file name
end

% Load CellBase
load(getpref('cellbase','fname'));
CELLIDLIST = CELLIDLIST(:);  %#ok<NODEF> % convert to column vector

% Write Excel file
TheMatrix_xls = formatforExcel(TheMatrix);   % change special numbers to strings
status1 = xlswrite(fname,CELLIDLIST,'sheet1','A1');   % export cellIDs
status2 = xlswrite(fname,TheMatrix_xls,'sheet1','B1');   % export TheMatrix
status = status1 & status2;   % status = true if both the cell IDs and TheMatrix has been written