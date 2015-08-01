function plot_raster(fighandle,time,binraster,trial_order,EventTimes,tlimits,partitions,partition_colors,labelx,labely);


% trial_order = valid_trials(ind);
% ShowEvents = {'OdorPokeIn','OdorValveOn','OdorPokeOut','WaterPokeIn','WaterPokeOut'};
% EventTimes = trialevents2relativetime(TE,g;TriggerEvent,ShowEvents);
% partitions{1} = part1;
% partitions{2} = part2;
% labelx=['Time-' g.TriggerEvent];
% labely='Firing rate';

trialNUM = length(trial_order):-1:1;

figure(fighandle)
clf;

hold on;
imagesc(time,trialNUM,1-binraster(trial_order,:));
colormap gray;

for iE = 1:size(EventTimes,1)
       plot(EventTimes(iE,trial_order),trialNUM, gcolor(iE,'.'))
end       
%plot(TE.OdorPokeIn(vind)   - refevent(vind),trialsREG,'y.')

axis([tlimits trialNUM(end) trialNUM(1)]);

plot_partitions(partitions{1},partition_colors{1},tlimits(1),10);
plot_partitions(partitions{2},partition_colors{2},tlimits(2),10);
     
xlabel(labelx);
ylabel(labely);
