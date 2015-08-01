function varargout = checkcb
%CHECKCB   Check main CellBase files.
%   CHECKCB checks if
%   - all TT files exist
%   - all Trialevents files exist
%   - all StimEvents files exist
%   - all STIMSPIKES files exist
%   - all EVENTSPIKES files exist
%
%   See also INITCB, LOADCB and PREALIGNSPIKES.

%   Edit log: BH 7/6/12

% Display message
disp('CHECKCB running');

% Check if CellBase exist
if isempty(getpref('cellbase')) 
     disp('CellBase does not exist.');
     return
end

% Compare CELLIDLIST to list of existing clusters 
cellids_listed = sort(listtag('cells'));    % defined by CellBase
cellids_exist = sort(findallcells);
if cellcmp(cellids_exist,cellids_listed)
    disp('Data found for all listed cellids.');
    CELLDATA_EXIST = 1;
else
    disp('Data and stored cellids don''t match.');
    CELLDATA_EXIST = 0;
end

% Check whether event files and prealigned files exist for all cells
PREALIGNED_STIMSPIKES_EXIST = 1;
PREALIGNED_EVENTSPIKES_EXIST = 1;
STIMEVENTS_EXIST = 1;
TRIALEVENTS_EXIST = 1;
for i = 1:length(cellids_listed)
    cellid = cellids_listed(i);
    cellstr = cellids_listed{i};
    fname1 = cellid2fnames(cellid,'StimEvent');
    fname2 = cellid2fnames(cellid,'STIMSPIKES');
    fname3 = cellid2fnames(cellid,'TrialEvent');
    fname4 = cellid2fnames(cellid,'EVENTSPIKES');
    
    if ~exist(fname1,'file');   % look for 'StimEvents'
        disp(['StimEvents data doesn''t exist for ' cellstr])
        STIMEVENTS_EXIST = 0;
    end
    if ~exist(fname2,'file');   % look for 'STIMSPIKES'
        disp(['Prealigned stimulation data doesn''t exist for ' cellstr])
        PREALIGNED_STIMSPIKES_EXIST = 0;
    end
    if ~exist(fname3,'file');   % look for 'TrialEvents'
        disp(['TrialEvents data doesn''t exist for ' cellstr])
        TRIALEVENTS_EXIST = 0;
    end
    if ~exist(fname4,'file');   % look for 'EVENTSPIKES'
        disp(['Prealigned behavioral data doesn''t exist for ' cellstr])
        PREALIGNED_EVENTSPIKES_EXIST = 0;
    end
end

% Display message if every events and prealigned files found
if TRIALEVENTS_EXIST
    disp('TrialEvents data exists for all cells.');
end
if PREALIGNED_EVENTSPIKES_EXIST
    disp('Prealigned behavior data exists for all cells.');
end
if STIMEVENTS_EXIST
    disp('StimEvents data exists for all cells.');
end
if PREALIGNED_STIMSPIKES_EXIST
    disp('Prealigned stimulation data exists for all cells.');
end

% Output = 1 if everything found
if nargout > 0
    varargout{1} = CELLDATA_EXIST & TRIALEVENTS_EXIST &PREALIGNED_EXIST;
end