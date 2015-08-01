1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     EXAMPLE SCRIPT FOR CELLBASE     %%%
%%%                                     %%%
%%%         AK, SPR 05,10/2010          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WORKFLOW FOR EACH SESSION (PREPROCESSING)
% use nlxcsc2mat2 to convert the CSCs into matlab use (do this on a Windows PC; ideally on the recording rig)
nlxcsc2mat2(pathname,'all');
% will convert all CSCs and Events to mat.
%% Cluster data from a session and copy entire session into database directory structure format.
% Base directory ) e.g. mPFCData
% ratname directory - e.g. d032 (Duda's animal32)
% recording session directories - yymmddx - e.g. 100407a (first recording session on 2010-04-07)
%% Initialize the database -- needs to be run once -- fast.
initcb
% Using spike times & event times prealigns and saves the data according to several criteria  -- run once, will take ~15 minutes.
allcells = listtag('cells'); % list all cellids
light_cellids=allcells([6 7 13 33 38 42 46 48 64 73 74 76 82 87 90 95 104 117 119 129 137 140 144 161 171 179 217 221 235]);
som_cellids=light_cellids(1:23);
pv_cellids=light_cellids(24:end);
%% Make sure everything is in the same timescale
match_timescales('Events')
match_timescales('spikes')
%% Run MakeTrialEvents and/or MakeStimEvents2 to convert events into a
%% trial based format
% pathname=uigetdir(getpref('cellbase','datapath'));
allcells=unique_session_cells;
problem_behav_cellid=[];
problem_stim_cellid=[];
allcells=listtag('cell');
for iC=1:length(allcells),
    cellid=allcells(iC);
    pathname=cellid2fnames(cellid,'Sess');
    try
        MakeTrialEvents2(pathname);
        script_behavior_analysis
        %pause
    catch
        problem_behav_cellid=[problem_behav_cellid cellid];
        disp(cellid);
    end
    try
        MakeStimEvents2(pathname,'BurstStartNttl',4)
    catch
        problem_stim_cellid=[problem_stim_cellid cellid];
    end
    try
        SE=load([pathname filesep 'StimEvents']);
        if isnan(SE.PulseOn),
            MakeStimEvents2(pathname,'BurstStartNttl',2)
        end
    catch
    end
end
% BurstStartNttl can be either 2 or 4 it was changed at some point but is
% now fixed.
%% prealignspikes with stimulation and behavioral events
% for iC=1:length(allcells),
tic;
allcells=listtag('cells');
allcells=allcells(97:111);
problem_behav_cellid=[];
problem_stim_cellid=[];
for iC=1:length(allcells),
    cellid=allcells(iC);
    disp(cellid)
    try
         %prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_backnforth3,'filetype','event') % filetype='event' for behavior protocol
         prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_NewEvents,'filetype','event') % filetype='event' for behavior protocol
    catch
        disp('Error in prealignSpikes.');
        problem_behav_cellid=[problem_behav_cellid cellid];
    end
%     try
%         prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_laserstim,'filetype','stim') % filetype='stim' for stimulation protocol
%     catch
%         disp(cellid)
%         problem_stim_cellid=[problem_stim_cellid cellid];
%     end
end
toc;
%% Extract cluster quality measures and waveforms
Calculate_CluSep
%% Extract waveform features for light activated and spontaneous waveforms
% saves it on database base directory
calc_light_spontaneous_waveform_features
%% View Event triggered rtaster PSTH
cellid={ 'd029_100119a_1.1' }  ;
TrigName='HomeZoneOut';
SEvent='TriggerZoneIn';

FNum=4;
BurstPSTH='off';
win=[-2 2];
parts='#SoundID';
parts='#PreviousRewardHistory';
parts='all';
dt=0.01;
sigma=0.02;
PSTHstd='on';
% % parts='#WaterValveDur';
ShEvent={{'TriggerZoneIn','RewardCue','TriggerZoneOut','Zone1FirstEntry','ReminderCue','Zone1FirstExit','WaterValveOn','WaterValveOff','HomeZoneIn','RewardZoneIn'}};
%ShEvent={{'RewardCue','TriggerZoneOut','ReminderCue','WaterValveOff','HomeZoneIn','RewardZoneIn','RewardZoneOut','HomeZoneOut','NextTriggerZoneIn'}};
% ShEvent={{'TriggerZoneIn','TriggerZoneOut'}};
ShEvColors=hsv(length(ShEvent{1}));
ShEvColors=mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
viewcell2b(cellid,'TriggerName',TrigName,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
'FigureNum',FNum,'eventtype','event','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
'EventMarkerWidth',0,'PlotZeroLine','on', BurstPSTH, 'BurstPSTH')
%% View Event triggered rtaster PSTH
%allcells=CELLIDLIST(69:end);
cellid={'d035_100503a_1.2'} ;
%allcells=allcells(223)
%allcells=som_cells(1:end);
TrigEvent='HomeZoneOut1';
SEvent='TriggerZoneIn';
FNum=7;
win=[-2 2];
parts='#SoundID'; 
parts='#PreviousRewardHistory'; 
parts='All';
dt=0.02; 
sigma=0.02;
PSTHstd='on';
% % parts='#WaterValveDur';
% ShEvent={{'TriggerZoneIn','RewardCue','TriggerZoneOut','Zone1FirstEntry','ReminderCue','Zone1FirstExit','WaterValveOn','WaterValveOff','HomeZoneIn','RewardZoneIn'}};
% ShEvent={{'RewardCue','TriggerZoneOut','ReminderCue','WaterValveOff','HomeZoneIn','RewardZoneIn','RewardZoneOut','HomeZoneOut','NextTriggerZoneIn'}};
ShEvent={{'TriggerZoneIn','RewardCue','ReminderCue','WaterValveOff','RewardZoneIn','RewardZoneOut','HomeZoneOut','NextTriggerZoneIn'}};
ShEvColors=hsv(length(ShEvent{1}));
ShEvColors=mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
% for iCell=1:length(allcells),
%     cellid=allcells(iCell);
    viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
    'FigureNum',FNum,'eventtype','event','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts)
% pause
% filename=['Light_stim_rasterpsth' '_' datestr(today) '_' cellid '.eps'];
% printfilename=[getpref('cellbase','datapath') filesep filename];
% end
%% View Light Triggered Raster PSTH
% eventtype stim
TrigEvent='BurstOn';
SEvent='BurstOff';
FNum=2;
win=[-0.005 0.01];
parts='all';
dt=0.001;
sigma=0.001;
PSTHstd='on';
ShEvent={{'PulseOn','PulseOff','BurstOff'}};
ShEvColors=hsv(length(ShEvent{1}));
ShEvColors=mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
% for iCell=1:length(allcells),
    cellid=allcells(1);
%     cellid=allcells(iCell);
%     try
        viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','on')
%         pause(1)
%     catch
%     end
% end

%% runanalysis
% Performs analysis and returns the values in a Matlab variable. It does
% not store it in the database
 mycellids=som_cells(end-10:end);
%mycellids=listtag('cells');
EpochName1='HomeZoneExit'; % immediately after HomeZoneExit
% EpochName2='RewardZoneExit';
EpochName2='PreHomeZoneOut'; % immediately before HomeZoenExit
nboot=1000;
eventtype='behav';
propvals=runanalysis(@calc_selectivity_epoch,'cellids',{mycellids},'property_names',{{'prop_name1'}},'arglist',{{EpochName1;EpochName2;nboot;eventtype}});
%% addanalysis
% Performs the analysis and adds the results to the database
% mycellids=pv_cellids;
mycellids=listtag('cells');
EpochName1='HomeZoneExit'; % immediately after HomeZoneExit
% EpochName2='RewardZoneExit';
EpochName2='PreHomeZoneOut'; % immediately before HomeZoenExit
nboot=1000;
eventtype='behav';
addanalysis(@calc_selectivity_epoch,'cellids',{mycellids},'property_names',{{'D_HZE_RZE','P_HZE_RZE','NT_HZE_RZE','RATE1_HZE_RZE','RATE2_HZE_RZE'}},'arglist',{{EpochName1;EpochName2;nboot;eventtype}});
%% viewpopselectivity

% Clean until now
%%
listtag('rat')      % lists all the rats in the dbase
listtag('session')  % lists all the sessions in the dbase

% Let's add an analysis 
% Runs an analysis on all units using the function 'meanrate'
addanalysis(@meanrate,'property_names',{{'MeanRate'}});

%Lists all the known properties, results of analyses
listtag('properties')   % or it's enough to type: listtag('prop')

%Lists the property 'Rate' for all units
rates = getvalue('MeanRate');

%% Lets plot a psth
TriggerEvent = 'PulseOn';
Window = [-0.2 2];
FigNum = 1;
viewpsth(mycell,'TriggerEvent',TriggerEvent,'window',Window,'FigureNum',FigNum)

%% Now let's add some more something more sophisticated analyses
EpochName = 'TriggerZoneHangOut';
nboot = 200;
addanalysis(@calc_selectivity_side,'property_names',{{'Dside_TriggerZoneHangOut';'Pside_TriggerZoneHangOut'}},'arglist',{{EpochName;nboot}});

%%
EpochName = 'AfterChoice';
nboot = 200;
addanalysis(@calc_selectivity_outcome,'property_names',{{'Dout_AC';'Pout_AC'}},'arglist',{{EpochName;nboot}});
%
EpochName1 = 'Baseline';
EpochName2 = 'AfterChoice';
nboot = 200;
addanalysis(@calc_selectivity_epoch,'property_names',{{'Dmod_AC';'Pmod_AC'}},'arglist',{{EpochName1;EpochName2;nboot}})
%
EpochName1 = 'OdorSampling';
EpochName2 = 'Baseline';
nboot = 0;
addanalysis(@calc_selectivity_epoch,'property_names',{{'Dmod_OSD';'Pmod_OSD'}},'arglist',{{EpochName1;EpochName2;nboot}})
%
%
listtag('anal')  % all the analyses ran
listtag('prop')  % all the property names
%
% Not the best function, but you can use it to delete the last analysis
% using a funhandle.
% delanalysis(@calc_selectivity_epoch);
%
mycell = allcells{1};
myvalues = listcell(mycell);   % lists all the property values associated with a cell

delcell(mycell);             %delete a cell
addcell(mycell);             %add a cell & execute all previous analyses


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do some trial selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is mostly internal to functions

% loadcb w/o arguments loads the database, with a cellid it loads the cell
% and you can further specify other files such as 'Events' associated with a cell
TE = loadcb(mycell,'Events');  % Loads the TrialEvents file associated with mycell

%You can select trial types of interest with this simple syntax
trialpos = selecttrial(TE,'OdorSampDur > 0.1 & OdorPokeValid & ChoiceDir == 2');  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's identify cells of interest:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% There might be better ways of defining search functions

%Locates cells based on rat / session / tetrode
myrat = findcell('rat','n1');
%Selects based on property names (very basic parser)
myselectOUT = selectcell('"Pout_AC" < 0.05 & "Dout_AC" < 0');       % significant & error outcome selective in AfterChoice period
myselectSIDE = selectcell('"Pside_MOV" < 0.05 & "Dside_MOV" < 0');   % significant & left selective in Movement period
mycellids = intersect(myrat,myselectOUT); 


%% Run an analysis without adding to the database or only on selected cells

RATES = runanalysis(@meanrate,'cellids',{mycellids}); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do some population analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%histogram of one property
viewpopselectivity(allcells,'Property','MeanRate');

%histogram with significance threshold
viewpopselectivity(allcells,'Property','Dout_AC','Significance','Pout_AC','Pthreshold',0.05);

%scatterplot of one property against another with significance coloring
viewpopselectivity(allcells,'Type','scatter','Property1','D_TriggerZoneHangOut_ReturnRunFirstPart','Significance1','P_TriggerZoneHangOut_ReturnRunFirstPart','Property2','D_TriggerZoneHangOut_ForwardRun','Significance2','P_TriggerZoneHangOut_ForwardRun');

%scatterplot as above but values are absolute (Transform - abs) --> now
%there is a nice correlates
viewpopselectivity(allcells,'Type','scatter','Property1','Dout_AC','Significance1','Pout_AC','Property2','Dmod_AC','Significance2','Pmod_AC','Transform1','abs','Transform2','abs');

%%
%plot some averages
som_prelimbic_ordered=som_prelimbic([  8 9 4  1 3  2 5 6 7 ]);
myselectOUT = som_prelimbic_ordered;

   %myselectOUT = som_cells(end-10:end);
  BurstPSTH='off';
  Partitions = {'Correct','Error'};
  Partitions = {'All'};
  % Partitions = {'PreviousRewardHistory'};
  % Partitions = {'SoundID'};
Sort= 'none';
Normalization ='zscore'; %mean, none etc.
TriggerEvent = 'HomeZoneOut1';
Window = [-5 5];
FigNum =13;
Overlay = 'on';
 viewpoppsth(myselectOUT,Partitions,'TriggerName',TriggerEvent,'Normalization',Normalization,...
    'window',Window,'FigureNum',FigNum,'Overlay',Overlay,'BurstPSTH',BurstPSTH);
viewpopraster(myselectOUT,Partitions,'TriggerName',TriggerEvent,'Normalization',Normalization,...
    'window',Window,'FigureNum',FigNum,'Overlay',Overlay,'BurstPSTH',BurstPSTH,'Sort',Sort);
% viewpoppsth(myselectOUT ,Partitions); % using all default values

%%
Partitions = {'ChoiceLeft','ChoiceRight'};  % Can be any variable in TrialEvents2 or 'MixtureDiff' or 'OdorRatio'
Normalization ='max';
TriggerEvent = 'OdorPokeOut';
Window = [-0.5 1];
FigNum = 2;
viewpoppsth(myselectSIDE,Partitions,'TriggerEvent',TriggerEvent,'Normalization',Normalization,'window',Window,'FigureNum',FigNum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Single unit analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
mycells = selectcell('"Pside_MOV" < 0.02 & "Dout_AC" < -0.1 & "Dmod_OSD" > 0.1');
mycell = mycells(1);

TriggerEvent = 'OdorPokeOut';
Window = [-0.2 2];
FigNum = 1;
viewpsth(mycell,'TriggerEvent',TriggerEvent,'window',Window,'FigureNum',FigNum)
%
TriggerEvent = 'WaterPokeIn';
SortEvent = 'WaterPokeOut';
Window = [-0.5 1.5];
FigNum = 1;
viewcell(mycell,'TriggerEvent',TriggerEvent,'SortEvent',SortEvent,'window',Window,'FigureNum',FigNum)
%
%
% a basic plot of the behavioral session for a cell
viewbehavior(mycells(4));

%%
%allcells=listtag('cells');

allcells=som_cells(26:end);
%allcells=allcells(223)
allcells={ 'd035_100501a_1.4'};
TrigName='HomeZoneOut1';
SEvent='TriggerZoneIn';
FNum=2;
NumPSTHPlots = 2;
LFPType = 'speed'
win=[-2 2];
parts='#SoundID';
 parts='#PreviousRewardHistory'; 
 parts='all';
dt=0.01; 
sigma=0.02;
PSTHstd='on';
% % parts='#WaterValveDur';
% ShEvent={{'TriggerZoneIn','RewardCue','TriggerZoneOut','Zone1FirstEntry','ReminderCue','Zone1FirstExit','WaterValveOn','WaterValveOff','HomeZoneIn','RewardZoneIn'}};
% ShEvent={{'RewardCue','TriggerZoneOut','ReminderCue','WaterValveOff','HomeZoneIn','RewardZoneIn','RewardZoneOut','HomeZoneOut','NextTriggerZoneIn'}};
ShEvent={{'TriggerZoneIn','RewardCue','ReminderCue','WaterValveOff','RewardZoneIn','RewardZoneOut','HomeZoneOut','NextTriggerZoneIn'}};
ShEvColors=hsv(length(ShEvent{1})); ShEvColors(1,:)=[0 1 1];
ShEvColors=mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
for iCell=1:length(allcells),
    cellid=allcells(iCell);
    viewcell2c_speed(cellid,'TriggerName',TrigName,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
'FigureNum',FNum,'eventtype','event','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
'EventMarkerWidth',0,'PlotZeroLine','on','NumPSTHPlots',2,'Num2Plot',20,'LFPType','speed')
pause
% filename=['Light_stim_rasterpsth' '_' datestr(today) '_' cellid '.eps'];
% printfilename=[getpref('cellbase','datapath') filesep filename];
end

%%
