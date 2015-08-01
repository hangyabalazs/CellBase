function initcb(cb_name)
%INITCB   Initialize CellBase.
%   INITCB specifies default directory and main database file for CellBase.
%   It also set CellBase preferences.
%
%   INITCB(NAME) creates new CellBase with the specified name. Multiple
%   instances of CellBase can be initialized with different names. Default
%   CellBase can be set using CHOOSECB.
%
%   See also LOADCB, CHOOSECB, DELETECB and ADDNEWCELLS.

%   Edit log: AK 3/04, 10/06; BH 3/18/11, 5/30/11, 4/26/12, 5/7/12, 8/20/13

% Check input arguments
error(nargchk(0,1,nargin))
if nargin < 1
    cb_name = '';
else
    cb_name = checknmcb(cb_name);
end

% Specify default variables
datapath = pwd;
fname = 'CellBase.mat';

% Set CellBase directory and main database file
d = uigetdir(datapath,'Select the root data directory');
if d ~= 0
   datapath = d;
   setpref('cellbase','datapath',datapath);
   [pathstr name ext] = fileparts(fname);
   fname = fullfile(datapath,[name ext]);
   f = uiputfile({'*.mat','MAT-files (*.mat)'},'Select the main database file',fname);
   if f ~=0
      fname = fullfile(datapath,f);
      if isempty(cb_name)
          cb_name = input('Give a name to your CellBase! ','s');
      end
      cb_name = checknmcb(cb_name);
      setpref('cellbase','name',cb_name);
      setpref('cellbase','fname',fname);
   else
      disp('INITCB canceled');
      return
   end
else
   disp('INITCB canceled');
   return
end
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
gp = getpref('cellbase');
if isfield(gp,'cellbases')
    gp = rmfield(gp,'cellbases');
end
cellbases{end+1} = gp;
setpref('cellbase','cellbases',cellbases)

% Open main database file
QuestionStr = sprintf('\n %s already exists. \n\n Do you want to delete it? \n',fname);
fid = fopen(fname);
if ( fid == -1) || strcmp(questdlg(QuestionStr,...
        'InitCB','Yes','No','No'),'Yes') %file doesn't exist or to be deleted !   
    
    % Save empty CellBase
    TheMatrix = [];
    ANALYSES  = [];
    CELLIDLIST = [];
    save(fname,'TheMatrix','ANALYSES','CELLIDLIST');
    
    % Convert TrialEvents & Sessions
    NUM_CELLS = addnewcells;
    funhandle = @cellid2vals;
    try
        addanalysis(funhandle,'property_names',{{'RatId';'DateNum';'Tetrode';'Unit'}});
        disp('Analysis added: cellid2vals.')
    catch
        disp('No analysis added.')
    end
%     if NUM_CELLS    % let's celebrate
        Data = imread('CellBase_icon.tif');
        welcomestr = sprintf('Your CellBase has been initialized with %d cells. \n Have fun!  \n\n\n\n',NUM_CELLS);
        msgbox(welcomestr,'Welcome to CellBase','custom',Data);
%     else
%         warndlg('Something went wrong.','InitCB');
%     end
else
    fclose(fid);
    welcomestr = sprintf('\n %s already exists. Delete it if you need a rebuild.\n',fname);
    msgbox(welcomestr,'Welcome to CellBase','warn');
end