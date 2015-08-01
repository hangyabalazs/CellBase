function viewpsth(cellid,varargin)
%
% VIEWPSTH
%
% -- not working flexibly yet


% trigger_event = 'WaterPokeIn'; 
% compare = {'Correct','Error','Left','Right'};
% window = [-0.5 1];
% dt = 0.01;
% sigma = 0.02;

if (nargin < 1)
	help viewpsth
	return
end

if  validcellid(cellid,{'list'}) ~= 1
    fprintf('%s is not valid.',cellid);
    return
    
end
% put arguments into a structure
% if ~isempty(varargin)
%     try, g = struct(varargin{:}); 
%     catch, error('Argument error in the {''param'', value} sequence'); end; 
% end;
% 
% % assign defaults if no arguments passed
% try, g.TriggerEvent;    catch, g.TriggerEvent = 'TrialStart'; end;
% try, g.SortEvent;       catch, g.SortEvent    = 'TrialStart';   end;
% try, g.Compare;         catch, g.Compare = {}; end;
% try, g.partition;       catch, g.partition     = '';  end;
% try, g.window;          catch, g.window        = [-0.0025 0.005];  end;
% try, g.dt;              catch, g.dt            = 0.001; end;
% try, g.sigma;           catch, g.sigma         = 0.002; end;
% try, g.plot;            catch, g.plot          = 'on'; end;
% try, g.FigureNum;       catch, g.FigureNum     = 1; end;

default_args={...
    'TriggerEvent',     'TrialStart';...
    'SortEvent',        'TrialStart';...
    'Compare',          {};...
    'filetype',         'event';... % 'event' 'stim' This determines what type of events we are prealigning to
    'partition',        '';...
    'window',           [-0.5 1];...
    'dt',               0.05;...
    'sigma',            0.1;...
    'plot',             'on';...
    'FigureNum',        1;...
    };
[g, error] = parse_args(default_args,varargin{:});

% % assign defaults if no arguments passed
% try, g.TriggerEvent;    catch, g.TriggerEvent = 'WaterPokeIn'; end;
% try, g.SortEvent;       catch, g.SortEvent    = 'WaterPokeOut';   end;
% try, g.Compare;         catch, g.Compare = {'Correct','Error','ChoiceLeft','ChoiceRight'}; end;
% try, g.partition;       catch, g.partition     = 'stimulus';  end;
% try, g.window;          catch, g.window        = [-0.25 1];  end;
% try, g.dt;              catch, g.dt            = 0.01; end;
% try, g.sigma;           catch, g.sigma         = 0.02; end;
% try, g.plot;            catch, g.plot          = 'on'; end;
% try, g.FigureNum;       catch, g.FigureNum     = 1; end;
%      
% test arguments for consistency
switch lower(g.plot)
    case { 'on', 'off' }, ;
    otherwise error('PLOT must be either on or off');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cellid='N49_050803_4.1';

TE = loadcb(cellid,'Trialevents');
ST = loadcb(cellid,'STIMES');

trigger_pos = findcellstr(ST.events(:,1),g.TriggerEvent);

if (trigger_pos == 0)
  error('Trigger variable not found');
end

%# valid_trials = find(~isnan(s.OdorPokeIn));  % where there was an OdorPokeIn at least
alltrials = 1:size(ST.event_stimes{1},2);

stimes  = ST.event_stimes{trigger_pos}(alltrials);;
windows = ST.event_windows{trigger_pos}(:,alltrials);

NUM_TRIALS = length(alltrials);

% add an extra margin to the windows
margin = g.sigma*3;

% time base array
time = g.window(1)-margin:g.dt:g.window(2)+margin;
    
%%% MAKE THE MAIN RASTER
binraster = stimes2binraster(stimes,time,g.dt);


% g.Compare = {'Correct','Error','ChoiceLeft','ChoiceRight','LeftCorrect','LeftError','RightCorrect','RightError','Stimulus'};
g.Compare={'StimType'};

for i = 1:size(g.Compare,2)  
  COMP(i,:) = int16(eval(['TE.' g.Compare{i}]));
end
%INC(1:length(TE.Correct)) = 1;     %include all
%INC(find(TE.Correct==0)) =  0; %exclude trials that were neither correct nor incorrect

% Include trials where the two major events are valid
% INC= TE.OdorPokeValid & TE.WaterPokeValid;
 INC=ones(1,size(COMP,2));
[psth, spsth, spsth_se] = binraster2psth(binraster,g.dt,g.sigma,COMP,INC);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=make3color([0 1 0],[0 0 1],[1 0 0]);
junk=m(1:11:end,:);
%colors6=mat2cell(junk,[1 1 1 1 1 1]);
colors6 = num2cell(junk',[1 3]);



labelx=['Time-' g.TriggerEvent];
labely='Trials';
fighandle = g.FigureNum;
tlimits = [g.window(1) g.window(2)];

figure(g.FigureNum)
clf;
col={'g','r',colors6{3},colors6{6},'c','m','k',colors6{:},'b','y',colors6{:}};
subplot(221)
   hold on;
   for i=1:2
     errorshade(time,spsth(i,:),spsth_se(i,:),col{i});
   end
   legend('Correct','Error',2);
   pos_disp=restrict(time,g.window(1),g.window(2));
   MX = (max(max(spsth(:,pos_disp))) +  max(max(spsth_se(:,pos_disp))))*1.1;
   alim=[g.window(1) g.window(2) 0 MX];
   axis(alim);
   xlabel(labelx); ylabel('Firing rate');
subplot(222)
   hold on;
   for i=3:4
     errorshade(time,spsth(i,:),spsth_se(i,:),col{i});
   end   
   legend('Left','Right',2);
   axis(alim);
   xlabel(labelx); ylabel('Firing rate');
subplot(223)
   hold on;
   for i=5:8
     errorshade(time,spsth(i,:),spsth_se(i,:),col{i});
   end   
   legend('LC','LE','RC','RE',2);
   axis(alim);
   xlabel(labelx); ylabel('Firing rate');
subplot(224)
   hold on;
   for i=9:size(spsth,1)
     %errorshade(time,spsth(i,:),spsth_se(i,:),col{i});
     plot(time,spsth(i,:),'color',col{i});
   end     
   legend(num2str([1:6]'),2)
   axis(alim)
   xlabel(labelx); ylabel('Firing rate');
