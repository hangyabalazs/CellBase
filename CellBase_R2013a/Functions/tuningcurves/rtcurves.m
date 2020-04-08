function rtcurves(cellids,varargin)
%RTCURVES   PSTHs for different reaction times.
%   RTCURVES(CELLIDS) calculates non-adaptive PSTHs and raster plots for
%   a set of cells aligned to 'Go' and 'No-go' tone/response onset. The
%   trials are partitioned based on reaction time quartiles. The PSTHs,
%   rasters and population average PSTHs are plotted and saved. Input
%   parameters:
%       CELLIDS - index set to CELLIDLIST (see CellBase documentation)
%
%   Default behavior can be modified by using a set of paramter-value pairs
%   as optional input parameters. The following parameters are implemented
%   (with default values):
%   	'align', 'tone' - align rtcurves to 'tone' or 'feedback'
%       'normalization', 'zscore' - use 'zscore' or 'max' normalizaztion;
%           default: Z-score normalization uses the mean across conditions
%           to obtain mean and SD (normalizing constants) to preserve 
%           differences of conditions.
%       'doraster', false - controls plotting of raster plots
%   	'issave', false - controls saving
%       'resdir', tempdir - results directory; uses system temporary folder
%           as default
%
%   See also NBPERFORMANCECURVES, NBTUNINGCURVES, ULTIMATE_PSTH and
%   NBRESPONSESORTER.

%   Edit log: BH 9/19/12, 2/13/14, 10/1/17
 
% Input argument check
cellids = cellids(:)';   % convert to row vector
prs = inputParser;
addRequired(prs,'cellids')  % list of cellIDs
addParameter(prs,'align','tone',@(s)ischar(s))   % controls alignment: 'tone' or 'response' or 'feedback'
addParameter(prs,'normalization','zscore',@(s)ischar(s))   % controls normalization: 'zscore' or 'max'
addParameter(prs,'doraster',false,@(s)islogical(s))   % controls displaying raster plots
addParameter(prs,'issave',false)   % controls saving
addParameter(prs,'resdir',tempdir,@isdir)   % results directory (default: system temp folder)
parse(prs,cellids,varargin{:})
g = prs.Results;

% Directories
if g.issave
    fs = filesep;
    resdir = [g.resdir g.align '_' g.normalization fs];   % results directory
    if ~isdir(resdir)
        mkdir(resdir)
    end
end

% PSTH
poppsth(cellids,g.align,'hit',g.normalization,g.doraster,g.issave,resdir);
poppsth(cellids,g.align,'fa',g.normalization,g.doraster,g.issave,resdir);

% -------------------------------------------------------------------------
function poppsth(I,align,trialtype,normalization,doraster,issave,resdir)

% Time window
wn = [-1.2 1.2];   % in seconds
dt = 0.001;
time = wn(1):dt:wn(2);   % time vector

% Call 'ultimate_psth'
% Preallocate
problem_ids = {};   % cellIDs for which PSTH failed
problem_msg = {};   % corresponding error messages
allpsth_orig = [];  % all PSTHs, for Hit and Miss
allspsth_orig = []; % all smoothed PSTHs, for Hit and Miss
NumCell = length(I);  % number of cells
for k = 1:NumCell
    cellid = I{k};
    disp(cellid)
    
    % Control PSTH and raster alignment
    switch align
        case 'tone'
            hitevent = 'StimulusOn'; %#ok<NASGU>
            faevent = 'StimulusOn'; %#ok<NASGU>
            sevent = 'StimulusOff';
        case 'response'
            hitevent = 'LeftWaterValveOn'; %#ok<NASGU>
            faevent = 'LeftPortIn'; %#ok<NASGU>
            sevent = 'StimulusOn';
        case 'feedback'
            hitevent = findAlignEvent_posfeedback_gonogo(cellid); %#ok<NASGU>
            faevent = findAlignEvent_negfeedback_gonogo(cellid); %#ok<NASGU>
            sevent = 'StimulusOn';
    end
    cevent = eval([trialtype 'event']);  % hitevent or faevent
    hitfilter = 'selectGoRT'; %#ok<NASGU>
    fafilter = 'selectNoGoRT'; %#ok<NASGU>
    cfilter = eval([trialtype 'filter']);   % hitfilter or fafilter
    
    try
        
        % Calcualte PSTH for Hits
        [psth1, spsth1, spsth_se1, ~, spt1] = ultimate_psth(cellid,'trial',cevent,wn,...
            'dt',dt,'display',true,'sigma',0.02,...
            'isadaptive',2,...
            'event_filter',cfilter,'filterinput',[0 0.25],...
            'maxtrialno',Inf,'baselinewin',[-0.5 0],'testwin',[0 0.5],...
            'relative_threshold',0.1);
        [psth2, spsth2, spsth_se2, ~, spt2] = ultimate_psth(cellid,'trial',cevent,wn,...
            'dt',dt,'display',true,'sigma',0.02,...
            'isadaptive',2,...
            'event_filter',cfilter,'filterinput',[0.25 0.5],...
            'maxtrialno',Inf,'baselinewin',[-0.5 0],'testwin',[0 0.5],...
            'relative_threshold',0.1);
        [psth3, spsth3, spsth_se3, ~, spt3] = ultimate_psth(cellid,'trial',cevent,wn,...
            'dt',dt,'display',true,'sigma',0.02,...
            'isadaptive',2,...
            'event_filter',cfilter,'filterinput',[0.5 0.75],...
            'maxtrialno',Inf,'baselinewin',[-0.5 0],'testwin',[0 0.5],...
            'relative_threshold',0.1);
        [psth4, spsth4, spsth_se4, ~, spt4] = ultimate_psth(cellid,'trial',cevent,wn,...
            'dt',dt,'display',true,'sigma',0.02,...
            'isadaptive',2,...
            'event_filter',cfilter,'filterinput',[0.75 1],...
            'maxtrialno',Inf,'baselinewin',[-0.5 0],'testwin',[0 0.5],...
            'relative_threshold',0.1);
        
        H = figure;
        psth = cat(1,psth1,psth2,psth3,psth4);
        spsth = cat(1,spsth1,spsth2,spsth3,spsth4);
        spsth_se = cat(1,spsth_se1,spsth_se2,spsth_se3,spsth_se4);
        tags = {'RT 1st quartile','RT 2nd quartile','RT 3rd quartile','RT 4th quartile'};
        
        % Concatenate data from different cells
        tagsn = extractnum(tags);  % get number from tag strings
        [st stia] = sort(tagsn);
        if isempty(allpsth_orig)
            linx = 1;  % next index
        else
            linx = size(allpsth_orig,2) + 1;
        end
        disp(size([spt1; spt2; spt3; spt4],1))  % display number of trials
        for tgs = 1:length(stia)
            allpsth_orig(tgs,linx,1:length(time)) = psth(stia(tgs),:); %#ok<AGROW>
            allspsth_orig(tgs,linx,1:length(time)) = spsth(stia(tgs),:); %#ok<AGROW>
        end
        
        % Plot PSTH
        NumPsth = size(spsth,1);
        clr = summer(NumPsth);   % colormap
        for pk = 1:NumPsth
            errorshade(time,spsth(pk,:),spsth_se(pk,:),'LineWidth',2,'LineColor',clr(pk,:),...
                'ShadeColor',clr(pk,:))
        end
        warning off
        legend(tags)
        warning backtrace
        maximize_figure(H)
        
        % Save PSTH figure
        if issave
            cellidt = regexprep(cellid,'\.','_');
            title([regexprep(cellidt,'_',' ') ' ' trialtype])
            fnm = [resdir cellidt '_PSTH_' trialtype '.fig'];
            saveas(H,fnm)
        end
        close(H)
        
        % Raster plot
        if doraster
            H = figure;
            viewcell2b(cellid,'TriggerName',cevent,'SortEvent',sevent,...
                'eventtype','behav','ShowEvents',{{sevent}},'PSTHstd','off',...
                'Partitions','#StimulusDuration&Hit','window',[-5 5])
            maximize_figure(H)
            
            % Save raster plot figure
            if issave
                fnm = [resdir cellidt '_raster_' trialtype '.fig'];
                saveas(H,fnm)
                fnm = [resdir cellidt '_raster_' trialtype '.jpg'];
                saveas(H,fnm)
            end
            close(H)
        end
        
    catch ME
        
        % Error handling
        problem_ids{end+1} = cellid;  %#ok<AGROW> % collect cellIDs resulting in error
        problem_msg{end+1} = ME.message;  %#ok<AGROW> % store corresponding error messages
        disp(['Something went wrong for cell ' cellid '.'])
        disp(ME.message)   % display error message
    end
end

% Normalization for averaging
[lent lenc lenp] = size(allspsth_orig);
switch normalization
    case 'zscore'
        allpsth = (allpsth_orig - repmat(mean(mean(allpsth_orig,1),3),[lent,1,lenp])) ...
            ./ repmat(std(mean(allpsth_orig,1),[],3),[lent,1,lenp]);   % standardize based on mean across conditions
        allspsth = (allspsth_orig - repmat(mean(mean(allspsth_orig,1),3),[lent,1,lenp])) ...
            ./ repmat(std(mean(allspsth_orig,1),[],3),[lent,1,lenp]);  % standardize based on mean across conditions
    case 'max'
        allpsth = (allpsth_orig - repmat(mean(mean(allpsth_orig,1),3),[lent,1,lenp])) ...
            ./ repmat(max(mean(allpsth_orig,1),[],3),[lent,1,lenp]);   % normalize based on maximum of mean across conditions
        allspsth = (allspsth_orig - repmat(mean(mean(allspsth_orig,1),3),[lent,1,lenp])) ...
            ./ repmat(max(mean(allspsth_orig,1),[],3),[lent,1,lenp]);  % normalize based on maximum of mean across conditions
end

% Normalized population average
allmn = squeeze(nanmean(allpsth,2));
allse = squeeze(nanstd(allpsth,[],2)) / sqrt(lenc);
allmns = squeeze(nanmean(allspsth,2));
allses = squeeze(nanstd(allspsth,[],2)) / sqrt(lenc);

% Plot population average
clr = summer(lent);
figure
hold on
for k = 1:lent
    errorshade(time,allmn(k,:),allse(k,:),'LineWidth',2,'LineColor',clr(k,:),...
        'ShadeColor',clr(k,:))
end
xlim([time(1) time(end)])
legend(tags)

% Save PSTH figure
H = gcf;
if issave
    fnm = [resdir 'meanPSTH_' trialtype '.fig'];
    saveas(H,fnm)
end
close(H)

% Plot population average
clr = summer(lent);
figure
hold on
for k = 1:lent
    errorshade(time,allmns(k,:),allses(k,:),'LineWidth',2,'LineColor',clr(k,:),...
        'ShadeColor',clr(k,:))
end
xlim([time(1) time(end)])
legend(tags)

% Save PSTH figure
H = gcf;
if issave
    fnm = [resdir 'meanSPSTH_' trialtype '.fig'];
    saveas(H,fnm)
end
close(H)

% keyboard

% -------------------------------------------------------------------------
function tagsn = extractnum(tags)

NumTags = length(tags);
tagsn = nan(1,NumTags);
for k = 1:NumTags
    intoken = tags{k};
    tkn = regexprep(intoken,'.*(\d).*','$1');   % extract number from the tag string
    tagsn(k) = str2double(tkn);
end