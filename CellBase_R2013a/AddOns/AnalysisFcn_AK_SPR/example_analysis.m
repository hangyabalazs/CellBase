%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     EXAMPLE SCRIPT FOR CELLBASE     %%%
%%%                                     %%%
%%%         AK 11/2006                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the dabase -- needs to be run once -- fast.
%
initcb
%
% Using spike times & event times prealigns and saves the data according
% several criteria  -- run once, will take ~15 minutes.
%
prealignSpikes('all');

allcells = listtag('cells');       % list all cellids

listtag('rat')      % lists all the rats in the dbase
listtag('session')  % lists all the sessions in the dbase

% Let's add an analysis 
% Runs an analysis on all units using the function 'meanrate'
addanalysis(@meanrate,'property_names',{{'MeanRate'}});

%Lists all the known properties, results of analyses
listtag('properties')   % or it's enough to type: listtag('prop')

%Lists the property 'Rate' for all units
rates = getvalue('MeanRate');

% Lets plot a psth
TriggerEvent = 'PulseOn';
Window = [-0.2 2];
FigNum = 1;
viewpsth(mycell,'TriggerEvent',TriggerEvent,'window',Window,'FigureNum',FigNum)

%Now let's add some more something more sophisticated analyses
EpochName = 'Movement';
nboot = 200;
addanalysis(@calc_selectivity_side,'property_names',{{'Dside_MOV';'Pside_MOV'}},'arglist',{{EpochName;nboot}});
%
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
viewpopselectivity(allcells,'Type','scatter','Property1','Dout_AC','Significance1','Pout_AC','Property2','Dmod_AC','Significance2','Pmod_AC');

%scatterplot as above but values are absolute (Transform - abs) --> now
%there is a nice correlates
viewpopselectivity(allcells,'Type','scatter','Property1','Dout_AC','Significance1','Pout_AC','Property2','Dmod_AC','Significance2','Pmod_AC','Transform1','abs','Transform2','abs');

%plot some averages
Partitions = {'Correct','Error'};
Normalization ='max'; %mean, none etc.
TriggerEvent = 'WaterPokeIn';
Window = [-0.5 1.5];
FigNum = 1;
viewpoppsth(myselectOUT,Partitions,'TriggerEvent',TriggerEvent,'Normalization',Normalization,'window',Window,'FigureNum',FigNum);
% viewpoppsth(mycellids,Partitions); % using all default values

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

