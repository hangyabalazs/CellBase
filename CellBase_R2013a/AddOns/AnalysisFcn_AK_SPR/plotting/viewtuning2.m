function viewtuning2(cellid,varargin) 

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
SP = loadcb(cellid,'EVENTSPIKES');

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
    USE_OdorPairs = 1;  
end

    
TUNING = unique(TE.OdorRatio);

for iT = 1:length(TUNING)
  
      pos = find(TE.OdorRatio==TUNING(iT));
      VALUE(iT) = mean(EpochRate(pos));
      
%     trialsC = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & Correct ',OdorPairID));
%     trialsE = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & Error ',OdorPairID));
%     trials_diff = findpos(TE.OdorRatio,[32 44 56 68]);
%     trialsCdiff = intersect(trials_diff,trialsC);
%     trialsEdiff = intersect(trials_diff,trialsE);
%     trialsL = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & ChoiceLeft ',OdorPairID));
%     trialsR = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & ChoiceRight ',OdorPairID));
%     
%     [stimC{OdorPairID}, FstimC{OdorPairID}, FstimC_se{OdorPairID}] = tuning_curve(EpochRate(trialsC),TE.OdorRatio(trialsC));
%     [stimE{OdorPairID}, FstimE{OdorPairID}, FstimE_se{OdorPairID}] = tuning_curve(EpochRate(trialsE),TE.OdorRatio(trialsE));
%     
    
%     [stim, FstimE, FstimE_se] = tuning_curve(EpochRate,TE.OdorRatio(trialsE));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure(g.FigureNum)
clf;


end