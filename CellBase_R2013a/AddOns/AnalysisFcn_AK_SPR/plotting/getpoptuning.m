function [NVALUES, PTAGS, NTAGS] = getpoptuning(cellids,Partitions,varargin)
%
%   VIEWPOPTUNING
%
%

% assign defaults if no arguments passed
default_args = { ...
        'EpochName'             'AfterChoice'; ...
        'Analysis'              {{'nanmean(data)';'nanstd(data)/sqrt(length(data))'}};...
        'OdorPairID'            1;  ...
        'Normalization'         'max'; ...
        'NormalizationType'     'processed'; ...
        'NormalizationTrials'   'all'; ...
        'ValidTrials'           ''; ...
        'SortOutput'            'value';...
        'Data2Use'              'EpochRate';...
        %         'TriggerEvent'     'WaterPokeIn'; ...
        %         'LastEvents',        '';...
    };

[g, error] = parse_args(default_args,varargin{:});
g.Partitions = Partitions;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumCells = length(cellids);

%figure out how many partitions we have this is 
%complicated by the fact that not all sessions will 
%have the same conditions and hence the same number 
%of partitions.
uniqsescells = unique_session_cells(cellids);

for iUC=1:length(uniqsescells)
    TE = loadcb(uniqsescells{iUC},'Events');
    [TRIALS, ATAGS{iUC}, RTAGS, NUM_TAGS{iUC}] = partition_trials(TE,g.Partitions);
    NumParts(iUC) = length(TRIALS);
end

[PTAGS, indPTAGS] = unique([ATAGS{:}]);
junk = [NUM_TAGS{:}];
NTAGS = junk(indPTAGS);
NumPartitions = length(PTAGS);

if strcmpi(g.SortOutput,'value')
    [j, indTAGS]=sort(NTAGS);
elseif strcmpi(g.SortOutput,'tags')
    [j, indTAGS]=sort(PTAGS);
else
    indTAGS = 1:NumPartitions; %no sorting which is really 'tags' thanks to 'unique'
end
PTAGS = PTAGS(indTAGS);
NTAGS = NTAGS(indTAGS);

NumVals = length(g.Analysis);


NVALUES = nan(NumCells,NumPartitions,NumVals);

%NVALUES = NumCell x NumPartitions x [EpochRate Analysis]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%d cells to process.\n',NumCells)
for iCELL = 1:NumCells
 
print_progress(iCELL,round(NumCells/100),5);

cellid = cellids{iCELL};
TE = loadcb(cellid,'Events');

switch lower(g.Data2Use)
    case 'epochrate'
        
        SP = loadcb(cellid,'EVENTSPIKES');
        epoch_pos = findcellstr(SP.epochs(:,1),g.EpochName);
        if (epoch_pos == 0)
            error('Trigger variable not found');
        end
        DATA = SP.epoch_rates{epoch_pos};
    otherwise
        DATA = eval(g.Data2Use);
end

%More flexibility needed


alltrials = 1:length(TE.TrialStart);
%alltrials = 1:size(SP.event_stimes{1},2);
%%stimes  = ST.event_stimes{trigger_pos}(alltrials);;
%%windows = ST.event_windows{trigger_pos}(:,alltrials);

[COMPTRIALS, TAGS, RTAGS, NUM_TAGS] = partition_trials(TE,g.Partitions);

%%% Could be put as an option
if g.OdorPairID == 0
    valid_trials = selecttrial(TE,sprintf('OdorConc == 100 & OdorPokeValid & WaterPokeValid %s',g.ValidTrials));
else
    valid_trials = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & OdorPokeValid & WaterPokeValid %s',g.OdorPairID,g.ValidTrials));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

values = get_tuning(DATA,g.Analysis,COMPTRIALS,valid_trials);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(g.NormalizationType,'processed')
    normdata = values;
    nrepdim = 1;
else
    normdata = DATA(:);
    nrepdim = size(values,2);
end

switch lower(g.Normalization)
    case 'none'
        NORMfactor = 1;
    case 'mean'
        NORMfactor =  nanmean(normdata);
    case 'median', 'medi'
        NORMfactor =  nanmedian(normdata);
    case 'max'
        NORMfactor =  max(normdata);
    case 'perc90'   
         NORMfactor = prctile(normdata,90);
%     case 'maxrate'
%         NORMfactor = max(EpochRate(union(trialsC,trialsE)));
    otherwise
        NORMfactor = 1;
end

%here is the key to figure out which parts were calculated...
posPARTS = match_list(TAGS,PTAGS);

NVALUES(iCELL,posPARTS,1:NumVals) = values./repmat(NORMfactor,size(values,1),nrepdim);

end
