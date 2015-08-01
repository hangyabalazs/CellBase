function   cellid = fname2cellid(fname)
%FNAME2CELLID    Convert filenames to cell IDs.
%   CELLID = FNAME2CELLID(FILENAME) converts valid filenames into cellids
%   or returns empty if it fails.
%
%   Valid filenames
%   (1) start with the default path or only include 'rat\session\unit.mat'
%   (2) unit 1 of tetrode 1 is called Sc1_1.mat 
%   (3) session name can only contain '.' or '_' characters but not both
%   (4) should be consistent across the database
%
%   See alse FINDALLCELLS.

%   ZFM additions:
%   - Uses the preference 'cell_pattern' to store a search term that selects
%   files corresponding to units. The default is 'Sc*.mat'. E.g.
%
%   setpref('cellbase','cell_pattern','Sc*.mat')
%
%   - Uses the preference 'group' to store an 'cut' subdirectory
%   under the session directory. E.g.
%
%   setpref('cellbase','group','Mike');

%   Edit log: BH 3/21/11

% Get cellbase preferences
cellbase_fname = getpref('cellbase','fname');
cellbase_path  = getpref('cellbase','datapath');
if ispref('cellbase','cell_pattern')
    cell_pattern = getpref('cellbase','cell_pattern');
else
    cell_pattern = 'Sc';
end

% Strip datpath from file
fn = char(strrep(fname,cellbase_path,''));
fs = filesep;

% Parse the filename (ratname\sessionname\analysisdir\)
[ratname,remain]  = strtok(fn,fs);
[session,remain]  = strtok(remain(2:end),fs);

% Added options for specifying an analysis directory below the session
% directory
if ispref('cellbase','group')
    % Strip analysis path
    [ad_junk,remain] = strtok(remain,fs);
end

% Extract tetrode unit
[tetrodeunit,ext] = strtok(remain(2:end),'.');
tu =  sscanf(tetrodeunit,[cell_pattern '%d_%d']);
pos_u = strfind(session,'_');
pos_p = strfind(session,'.');

% Control output
if ~strcmp(ext,'.mat')
    %disp('FNAM2CELLID: Not a matlab data file.');
    cellid = 0;
    return
end

if isempty(ratname) || isempty(session) || isempty(tu)
   % disp('FNAM2CELLID: Filename could not be parsed correctly.') 
   cellid = 0;
   return
elseif ~isempty(pos_u) && ~isempty(pos_p)
    disp('FNAM2CELLID: Cannot use both periods and underscore in session names.');
    cellid = 0;
    return    
elseif ~isempty(pos_u)   % there were underscores in the sessions
    session = strrep(session,'_','.');   % replace them with .'s
    setpref('cellbase','session_separator','_');  % note this as a preference
elseif ~isempty(pos_p)
    setpref('cellbase','session_separator','.');
else
    % there is no separator, which is OK
    setpref('cellbase','session_separator','');
end

cellid = sprintf('%s_%s_%d.%d',ratname,session,tu(1),tu(2));
%disp(cellid)