function viewcell(cellid,varargin)
%
% VIEWCELL
%
%   -- not working very flexibly yet


if (nargin < 1)
	help viewcell
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
try, g.TriggerEvent;    catch, g.TriggerEvent = 'WaterPokeIn'; end;
try, g.SortEvent;       catch, g.SortEvent    = 'WaterPokeOut';   end;
try, g.partition;       catch, g.partition     = 'stimulus';  end;
try, g.window;          catch, g.window        = [-0.25 1];  end;
try, g.dt;              catch, g.dt            = 0.01; end;
try, g.sigma;           catch, g.sigma         = 0.02; end;
try, g.plot;            catch, g.plot          = 'on'; end;
try, g.FigureNum;       catch, g.FigureNum     = 1; end;
     
% test arguments for consistency
switch lower(g.plot)
    case { 'on', 'off' }, ;
    otherwise error('PLOT must be either on or off');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TE = loadcb(cellid,'Events');
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

% From this one can easily derive many useful quantities:

% % 
% % 
% % % DEAL WITH VARIABLE EPOCH WINDOWS 
% % min_trials = 2;
% % for iTIME=1:length(time)
% %     psth_trials(iTIME) = sum(windows(1,:) <= time(iTIME) & windows(2,:) > time(iTIME));
% %     % might want to watch out here
% %     if psth_trials(iTIME) < min_trials
% %         psth_trials(iTIME) = NaN;
% %     end
% % end

compare = {'Correct','Error','ChoiceLeft','ChoiceRight','LeftCorrect','LeftError','RightCorrect','RightError','Stimulus'};

for i = 1:size(compare,2)  
  COMP(i,:) = int16(eval(['TE.' compare{i}]));
end
%INC(1:length(TE.Correct)) = 1;     %include all
%INC(find(TE.Correct==0)) =  0; %exclude trials that were neither correct nor incorrect

% Include trials where the two major events are valid
INC = TE.OdorPokeValid & TE.WaterPokeValid;
 
[psth, spsth, spsth_se] = binraster2psth(binraster,g.dt,g.sigma,COMP,INC);

valid_trials = find(INC); 
partition2 = 'Outcome';
[ind, part1, part2] = trialsort(TE,g.partition,partition2,g.TriggerEvent,g.SortEvent,valid_trials);
trial_order = valid_trials(ind);
partitions{1} = part1;
partitions{2} = part2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=make3color([0 1 0],[0 0 1],[1 0 0]);
junk=m(1:11:end,:);
%colors6=mat2cell(junk,[1 1 1 1 1 1]);
colors6 = num2cell(junk',[1 3]);

switch lower(g.partition)
    case 'stimulus'
               colors1 = colors6;
    case 'direction'
               colors1 = {'b','r'};
end

colors2 = {'g','r'};

partition_colors{1} = colors1;
partition_colors{2} = colors2;


ShowEvents = {'OdorPokeIn','OdorValveOn','OdorPokeOut','WaterPokeIn','WaterPokeOut'};
EventTimes = trialevents2relativetime(TE,g.TriggerEvent,ShowEvents);

labelx=['Time-' g.TriggerEvent];
labely='Trials';
fighandle = g.FigureNum;
tlimits = [g.window(1) g.window(2)];

plot_raster(fighandle,time,binraster,trial_order,EventTimes,tlimits,partitions,partition_colors,labelx,labely);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   switch lower(compare{1})
%        case 'correct'
%            col={'g','r'};
%        case  'left'
%            col=colors6([3 6]);
%        case 'odorvalveid'
%            col = colors6;
%    end
figure(g.FigureNum+1)
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
   
figure(g.FigureNum+2)
clf;
%col={'g','r',colors6{3},colors6{6},'c','c--','m','m--',colors6{:}};
subplot(211)
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
subplot(212)
   hold on;
   for i=9:size(spsth,1)
     %errorshade(time,spsth(i,:),spsth_se(i,:),col{i});
     plot(time,spsth(i,:),'color',col{i});
   end     
   legend(num2str([1:6]'),2)
   axis(alim)
   xlabel(labelx); ylabel('Firing rate');