
function ClusterChecker (clname,sessionpath)

%Output a pdf 

if nargin ==0;
    [clname,sessionpath] = uigetfile(fullfile('.clusters'), 'Select a .clusters file'); 
end  
    if sessionpath(end)=='\'; [sessionpath]=fileparts(sessionpath);end
	clpath = fullfile(sessionpath, clname); %.clusters path name
    tname = strtok(clname, '.clusters');    %TT* file name
    [ratpath sessionname]=fileparts(sessionpath);
    [parentdir ratname]=fileparts(ratpath);
    
%% DefineParameters
    global MClust_NeuralLoadingFunction
    global MClust_FDdn MClust_FDext MClust_ChannelValidity MCLUST_DEBUG
    global MClust_TTfn MClust_TTdn MClust_TText MClust_TTData
    global MClust_Clusters MClust_FeatureTimestamps MClust_FeatureSources
    global MClust_ClusterSeparationFeatures
    
    %Variables
    pdfname = 'clusters.pdf';
    [t1 t2] = strtok(fliplr(which('MClust')),'\');
    MClust_Directory = fliplr(t2);
%     MClust_FDdn = [sessionpath '\FD\'];
    MClust_FDdn = [sessionpath];
    featureList = {'Amplitude'; 'Energy'; 'Time'; 'WavePC1'}; %% {'Amplitude'; 'Energy'; 'Time'; 'WavePC1'}
    MClust_ClusterSeparationFeatures = {'Energy' ;'WavePC1'}; %% {'Energy' ;'WavePC1'}
    MClust_ChannelValidity = [1 1 1 1];
    MClust_TTdn=sessionpath;
    
    %Constant
%     MClust_NeuralLoadingFunction=char([MClust_Directory 'LoadingEngines\LoadTT_NeuralynxNT']);
    MClust_NeuralLoadingFunction=char([MClust_Directory 'LoadingEngines\LoadTT_Intan']);
%     MClust_TText = '.ntt';
    MClust_TText = '.mat';
    MClust_FDext = '.fd';
    MClust_TTData = [];      % data from tt file
    MClust_FeatureSources = {}; % <filenames, number pairs> for finding features in fd files
    MClust_FeatureTimestamps=[];
    MCLUST_DEBUG=1;
    MClust_Clusters=[];
    MClust_TTfn=[];


%% LoadClusterData
	temp = load(clpath,'-mat');
    MClust_Clusters=temp.MClust_Clusters;
	nClusters = length(MClust_Clusters);
    
%% LoadFeatureData
    featureFiles =  sortcell(FindFiles('feature_*.m', ...
    'StartingDirectory', fullfile(MClust_Directory, 'Features') ,'CheckSubdirs', 0));

%% MakeFigures

MClust_TTfn=tname;
for iClust=1:nClusters;
CalculateFeatures(MClust_TTfn, featureList)

% Fix channel validity
tt = load([fullfile(MClust_FDdn, [MClust_TTfn '_' featureList{1}]) MClust_FDext],'-mat');
fdd = tt.FeatureData;
tmx = max(fdd);
MClust_ChannelValidity(tmx==0) = 0;
    


CO_01_CheckCluster(iClust)
% text(30,2000,[ratname '_' sessionname] ,'HorizontalAlignment', 'center','VerticalAlignment', 'top', 'fontsize', 20)
fstamp ([ratname '_' sessionname]);
set(gcf, 'Color', 'w');

export_fig(gcf,'-append',[ratname '_' sessionname '_' pdfname], '-zbuffer'); % just filename (cd) or specify fullpath

clf(gcf); close; %clear figure for each time after saving to release memory
end

