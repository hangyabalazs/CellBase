function x = stimes2selectivitytime(stimes,time,dt,ev_windows,window_margin);





g.TriggerEvent = 'WaterPokeIn';
cellid='N49_050803_4.1';

TE = loadcb(cellid,'Events');
SP = loadcb(cellid,'EVENTSPIKES');

trigger_pos = findcellstr(SP.events(:,1),g.TriggerEvent);

if (trigger_pos == 0)
  error('Trigger variable not found');
end

%# valid_trials = find(~isnan(s.OdorPokeIn));  % where there was an OdorPokeIn at least
alltrials = 1:size(SP.event_stimes{1},2);

stimes  = SP.event_stimes{trigger_pos}(alltrials);;

%-----------------------------
win
dwin

 | |   | | ||| ||   | |  |  | ||| | | 
 
 cumsum(x:y)
 