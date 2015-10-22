function mountcb(cb_name)
%MOUNTCB   Mount CellBase.
%   MOUNTCB mounts an existing CellBase in Matlab. The new CellBase is
%   added to the Matlab CellBase preferences.
%
%   MOUNTCB(NAME) mounts new CellBase with the specified name. Multiple
%   instances of CellBase can be initialized with different names. Default
%   CellBase can be set using CHOOSECB.
%
%   See also INITCB, LOADCB, CHOOSECB and DELETECB.

%   Edit log: TL 10/21/2015

% Check input arguments
narginchk(0,1)
if nargin < 1
    cb_name = input('Give a name to your CellBase! ','s');
end
cb_name = checknmcb(cb_name);  % force unique name

% Locate CellBase
[fname, datapath] = uigetfile({'*.mat','MAT-files (*.mat)'}, 'Locate CellBase database file');
if isequal(fname,0)
   disp('MOUNTCB canceled')
   return
else
   disp(['You selected ', fullfile(datapath, fname)])
end

% Set CellBase directory and main database file
fname = fullfile(datapath,fname);   % build the filename of the mat file - main database file
setpref('cellbase','datapath',datapath);   % sets CellBase preferences
setpref('cellbase','name',cb_name);
setpref('cellbase','fname',fname);
disp(['Data path is ' datapath]);
disp(['File name is ' fname])
disp(['CellBase name is ' cb_name])

% Timestamp conversion
tsc = questdlg('I''m going convert timestamps to seconds in','Timestamp conversion',...
    'CellBase','MClust','CellBase');
switch tsc
    case 'CellBase'
        timefactor = 1e-4;  % TT*.mat files will reflect the timestamps of the Ntt files;
                            % loadcb will convert them to seconds
                            % timefactor is 1 millisecond
    case 'MClust'
        timefactor = 1;     % timestamps in TT*.mat files will already be converted to seconds
end

% Set other preferences (persistent and maintain their values between MATLAB sessions)
setpref('cellbase','session_filename','TrialEvents.mat');
setpref('cellbase','TrialEvents_filename','TrialEvents.mat');
setpref('cellbase','StimEvents_filename','StimEvents.mat');
setpref('cellbase','cell_pattern','TT');
setpref('cellbase','filesep',filesep);
setpref('cellbase','timefactor',timefactor);

% Store cellbases to allow multiple instances
if ispref('cellbase','cellbases')
    cellbases = getpref('cellbase','cellbases');
else
    cellbases = {};
end
gp = getpref('cellbase');  % returns the preferences as a structure 
if isfield(gp,'cellbases')
    gp = rmfield(gp,'cellbases');
end
cellbases{end+1} = gp;  % adds the new CellBase to 'cellbases' 
setpref('cellbase','cellbases',cellbases)

% Feedback
welcomestr = sprintf('New CellBase named %s was successfully mounted.\n',cb_name);
msgbox(welcomestr,'Welcome to CellBase');