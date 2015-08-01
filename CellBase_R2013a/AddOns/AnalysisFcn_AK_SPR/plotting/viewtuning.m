function viewtuning(cellid,varargin) 

% 
% epoch_name = 'AfterChoice';
% nbins = 5;
% min_points = 10;
% PLOT = 4;
% %function viewcell(ev_stimes,ev_windows,events,s,trigger_event,sort_event,partition1,window,dt,sigma,PLOT)
% function viewcell(cellid,varargin)

if (nargin < 1)
	help viewtuning
	return
end

if  validcellid(cellid,{'list'}) ~= 1
    fprintf('%s is not valid.',cellid);
    return
    
end
% put arguments into a structure
if ~isempty(varargin)
    try, g = struct(varargin{:}); 
    catch, error('Argument error in the {''param'', value} sequence'); end; 
end;

% assign defaults if no arguments passed
try, g.EpochName;       catch, g.EpochName     = 'AfterChoice'; end;
try, g.Bins;            catch, g.Bins          = 5;   end;
try, g.MinPoints;       catch, g.MinPoints     = 10;  end;
try, g.plot;            catch, g.plot          = 'on'; end;
try, g.FigureNum;       catch, g.FigureNum     = 1; end;
try, g.OdorPairID;      catch, g.OdorPairID    = 0; end;
try, g.Pthreshold;      catch, g.Pthreshold    = 0.05; end;
     
% test arguments for consistency
switch lower(g.plot)
    case { 'on', 'off' }, ;
    otherwise error('PLOT must be either on or off');
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TE = loadcb(cellid,'Events');
SP = loadcb(cellid,'SPIKES');

% %# valid_trials = find(~isnan(s.OdorPokeIn));  % where there was an OdorPokeIn at least
% alltrials = 1:size(SP.event_stimes{1},2);
% 
% stimes  = SP.event_stimes{trigger_pos}(alltrials);;
% windows = SP.event_windows{trigger_pos}(:,alltrials);

epoch_pos = findcellstr(SP.epochs(:,1),g.EpochName);
if (epoch_pos == 0)
  error('Epoch name not found');
end

EpochRate = SP.epoch_rates{epoch_pos};

%valid_trials = find(TE.OdorPokeValid & TE.WaterPokeValid);


if g.OdorPairID == 0
    USE_OdorPairs = unique(TE.OdorPairID);
else
    USE_OdorPairs = unique(TE.OdorPairID);  
end

    
for OdorPairID = USE_OdorPairs
    Nmixes = length(unique(TE.OdorRatio(find(TE.OdorPairID == OdorPairID))));
    if Nmixes < 4
        USE_OdorPairs(OdorPairID) = []; % delete b/c no mixtures
        continue
    end
    
    trialsC = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & Correct ',OdorPairID));
    trialsE = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & Error ',OdorPairID));
    trials_diff = findpos(TE.OdorRatio,[32 44 56 68]);
    trialsCdiff = intersect(trials_diff,trialsC);
    trialsEdiff = intersect(trials_diff,trialsE);
    trialsL = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & ChoiceLeft ',OdorPairID));
    trialsR = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & ChoiceRight ',OdorPairID));
    
    [stimC{OdorPairID}, FstimC{OdorPairID}, FstimC_se{OdorPairID}] = tuning_curve(EpochRate(trialsC),TE.OdorRatio(trialsC));
    [stimE{OdorPairID}, FstimE{OdorPairID}, FstimE_se{OdorPairID}] = tuning_curve(EpochRate(trialsE),TE.OdorRatio(trialsE));
    
    
    [cafO{OdorPairID}, binsO{OdorPairID}, nt] = histratio(EpochRate(trialsC),EpochRate(trialsE),g.Bins,'fixed_bins','minpoints',g.MinPoints); 
    cafO_se{OdorPairID} = sqrt(cafO{OdorPairID}-cafO{OdorPairID}.^2) ./ sqrt(max(1,nt-1));  %binomiaml errorbar
    [DO(OdorPairID),PO(OdorPairID)] = rocarea(EpochRate(trialsC),EpochRate(trialsE),'boot',100,'scale');
    %
    [cafOdiff{OdorPairID}, binsOdiff{OdorPairID}, nt] = histratio(EpochRate(trialsCdiff),EpochRate(trialsEdiff),g.Bins,'fixed_bins','minpoints',g.MinPoints); 
    cafOdiff_se{OdorPairID} = sqrt(cafOdiff{OdorPairID}-cafOdiff{OdorPairID}.^2) ./ sqrt(max(1,nt-1));  %binomiaml errorbar
    [DOdiff(OdorPairID),POdiff(OdorPairID)] = rocarea(EpochRate(trialsC),EpochRate(trialsE),'boot',100,'scale');
    %
    [cafDIR{OdorPairID}, binsDIR{OdorPairID}, nt] = histratio(EpochRate(trialsL),EpochRate(trialsR),g.Bins,'fixed_bins','minpoints',g.MinPoints); 
    cafDIR_se{OdorPairID} = sqrt(cafDIR{OdorPairID}-cafDIR{OdorPairID}.^2) ./ sqrt(max(1,nt-1));  %binomiaml errorbar
    [DDIR(OdorPairID),PDIR(OdorPairID)] = rocarea(EpochRate(trialsL),EpochRate(trialsR),'boot',100,'scale');
    %
end

NumOdorPairs = length(USE_OdorPairs);
if NumOdorPairs < 1
    return
end
STAR = ' *';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure(g.FigureNum)
clf;

for OdorPairID = USE_OdorPairs

PlotNum = 3*(find(USE_OdorPairs == OdorPairID)-1);  % USE_OdorPairs is could be [1 3], so we need make sure the PlotNum is OK
    
subplot(NumOdorPairs,3,1+PlotNum)
hold on;
eC= errorbar(stimC{OdorPairID},FstimC{OdorPairID},FstimC_se{OdorPairID},'g');
eE= errorbar(stimE{OdorPairID},FstimE{OdorPairID},FstimE_se{OdorPairID},'r');
x=[FstimC{OdorPairID} FstimE{OdorPairID}];
axis([min(stimC{OdorPairID})-10 max(stimC{OdorPairID})+10 min(x)*0.8 max(x)*1.2]);
set([eE(2) eC(2)],'LineWidth',2);
set(gca,'XTick',stimC{OdorPairID})
xlabel('Stimulus'); ylabel('Firing rate');


subplot(NumOdorPairs,3,2+PlotNum)
hold on;
errorshade(binsO{OdorPairID},cafO{OdorPairID},cafO_se{OdorPairID},[0.8 0.8 0])
plot(binsO{OdorPairID},cafO{OdorPairID},'k*')
if ~isnan(binsOdiff{OdorPairID})
    errorshade(binsOdiff{OdorPairID},cafOdiff{OdorPairID},cafOdiff_se{OdorPairID},'b'); plot(binsOdiff{OdorPairID},cafOdiff{OdorPairID},'b*');
end
xlabel('Firing rate (Hz)'); ylabel('Accuracy');
title(['Outcome selectivity: ' num2str(DO(OdorPairID),2) STAR((PO(OdorPairID)<g.Pthreshold)+1)]);


subplot(NumOdorPairs,3,3+PlotNum)
hold on;
errorshade(binsDIR{OdorPairID},cafDIR{OdorPairID},cafDIR_se{OdorPairID},[0.8 0 0.8])
plot(binsDIR{OdorPairID},cafDIR{OdorPairID},'r*')
xlabel('Firing rate (Hz)'); ylabel('Direction');
title(['Direction selectivity: ' num2str(DDIR(OdorPairID),2) STAR((PDIR(OdorPairID)<g.Pthreshold)+1)]);


end