function allcellids = findallcells(varargin)
%FINDALLCELLS   Find all clustered cells on the CellBase path.
%   IDS = FINDALLCELLS finds all clustered cells in the default cellbase
%   directory and subdirectories.
%
%   IDS = FINDALLCELLS(DR) returns cell IDs from DR.
%
%   See also FINDCELLPOS.

%   Edit log: AK  3/04; ZFM 7/05; BH 3/21/11, 4/23/13

% Get cellbase preferences
cellbase_datapath = getpref('cellbase','datapath');

% Input argument check
if nargin > 0
    ratdir = varargin{1};
    if ~exist(fullfile(cellbase_datapath,ratdir),'dir')
        disp('Directory doesn''t exist.')
        return
    end
else
    ratdir = listdir(cellbase_datapath);
end
if ischar(ratdir)
    ratdir = {ratdir};
end
if isempty(ratdir)
    disp('FINDALLCELLS: No cells found.');
    return
end

% Get cell IDs
k = 1;  % counter for cell IDs
allcellids = {};
for rdir = 1:length(ratdir)   % animal loop
    ratpath = fullfile(cellbase_datapath,char(ratdir(rdir)));
    sessiondir = [{''} listdir(char(ratpath))];
    for sdir = 1:length(sessiondir)   % session loop
        
        % Use cell_pattern property
        fullsessiondir = fullfile(ratpath,char(sessiondir(sdir)));
        if ispref('cellbase','cell_pattern')
            cell_pattern = getpref('cellbase','cell_pattern');
            
            % Check for cells
            cellfiles = listfiles(fullsessiondir,[cell_pattern '(\d)_(\d).mat']);
        else
            cellfiles = listfiles(fullsessiondir,'TT(\d)_(\d).mat');
        end
        
        % Convert filenames to cell IDs
        for fnum = 1:length(cellfiles)   % filename loop
            fname = fullfile(fullsessiondir,char(cellfiles(fnum)));
            cellid = fname2cellid(fname);
            if cellid ~= 0
                allcellids{k} = cellid;
                k = k+1;
            end
        end   % end of filename loop
    end   % end of session loop
end   % end of animal loop