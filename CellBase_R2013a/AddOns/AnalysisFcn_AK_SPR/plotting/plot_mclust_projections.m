function plot_mclust_projections(cellid)
% PLOT_MCLUST_PROJECTIONS scatter plot of optical tagging cluster
% projection. Plots a scatter plot of energy values (Energy 2 versus Energy
% 4) calculated from mclust. All neurons recorded on the tetrode are shown
% and light evoked spikes are overlayed on optically tagged cluster in
% blue. Only 10% of noise spikes and 50% of clustered spikes are plotted.
% Light evoked spikes are selected in a window of 2 ms from pulse onset.
% Only the first light evoked spike is selected for plotting. To account
% for drift, only spikes during the stimulation protocol are selected for
% plotting.
% SPR 2011-12-28

if nargin < 1,
    help plot_mclust_projections
    cellid = 'd033_100417a_5.2';
end

% Load stim events to get Pulse onset times.
ST = loadcb(cellid,'Stimevents');
pon = ST.PulseOn;
% pon = ST.BurstOn(~isnan(ST.BurstOn));
% pon = setdiff(ST.PulseOn,ST.BurstOn(~isnan(ST.BurstOn)));

% Cells on same tetrode including cellid.
[numcell,tetpartners] = tetrodepairs(cellid);
tag_cell = strmatch(cellid,tetpartners);

% decide colormap for all neurons.
mk_color = jet(numcell) + 0.2*ones(numcell,3);
mk_color(mk_color > 1) = 1;
mk_color = [0.1 0.1 0.1 ;mk_color];
mk_color = repmat(linspace(0.1,0.15,numcell+1)',1,3);

Nttfn = cellid2fnames(cellid,'Ntt');
% Load spikes from Ntt file.
all_spikes = LoadTT_NeuralynxNT(Nttfn);
all_spikes = all_spikes*1e-4;
nspk = length(all_spikes);
% consider spikes only within the stimulation protocol to account for
% drift.
[junkTS,wv]=LoadTT_NeuralynxNT(Nttfn,all_spikes*1e4,1); % Specified spikes

val_spk_i = [find(all_spikes >= pon(1),1,'first') find(all_spikes >= pon(end),1,'first')];

% Load feature data for tetrode.
[r,s,t,u] = cellid2tags(cellid);
prop = 'Energy.fd';
% prop = 'Peak.fd';

propfn = [getpref('cellbase','cell_pattern') num2str(t) '_' prop];
propfn_path = [cellid2fnames(cellid,'sess') filesep propfn];
wf_prop = load(propfn_path,'-mat');
wf_prop1 = wf_prop.FeatureData(:,1);
wf_prop2 = wf_prop.FeatureData(:,2);
wf_prop3 = wf_prop.FeatureData(:,3);
wf_prop4 = wf_prop.FeatureData(:,4);

% Spikes from each cell have an index. Spikes from noise have index 0.
cell_i = zeros(nspk,1);
for iCell = 1:numcell,
    % Load spike times.
    spk = loadcb(tetpartners(iCell),'Spikes');
    % get indices for the cell.
    [junk,junk2,cell_inx] = intersect(spk,all_spikes);
    cell_i(cell_inx) = iCell;
end

% Figure with 1 projection (2 Vs 4)
f1h = figure;
setmyfigure(gcf)
ax1 = axes;

% Figure with all projections plotted.
f2h = figure;
sub_h = set_subplots(2,3);

% Plot each cluster.
co = 0;
for iCell = 0:numcell,
    co = co+1;
    cell_inx = find(cell_i == iCell);
    cell_inx = cell_inx(cell_inx > val_spk_i(1) & cell_inx < val_spk_i(2));
    rnd_inx = randperm(length(cell_inx));
    if iCell == 0,
        rnd_inx = rnd_inx(1:ceil(0.1*length(cell_inx)));
    else
        rnd_inx = rnd_inx(1:ceil(0.5*length(cell_inx)));
    end
    cell_inx = cell_inx(rnd_inx);
    axes(ax1);
    if iCell == tag_cell,
        h = plot(wf_prop2(cell_inx),wf_prop4(cell_inx),'MarkerFaceColor',mk_color(co,:),'Marker','o','MarkerSize',4);
    elseif iCell == 0,
        h = plot(wf_prop2(cell_inx),wf_prop4(cell_inx),'MarkerFaceColor',mk_color(co,:),'Marker','o','MarkerSize',4);
    else
        h = plot(wf_prop2(cell_inx),wf_prop4(cell_inx),'MarkerFaceColor',mk_color(co,:),'Marker','o','MarkerSize',4);
    end
    hold on
    set(h,'LineStyle','none','MarkerEdgeColor','none')
    
    % Plot all Energy plots.
    figure(f2h)
    iS = 1;
    subplot(sub_h(iS));iS = iS+1;
    h = plot(wf_prop1(cell_inx),wf_prop2(cell_inx),'MarkerFaceColor',mk_color(co,:),'Marker','o','MarkerSize',4);
    hold on
    set(h,'LineStyle','none','MarkerEdgeColor','none')
    subplot(sub_h(iS));iS = iS+1;
    h = plot(wf_prop1(cell_inx),wf_prop3(cell_inx),'MarkerFaceColor',mk_color(co,:),'Marker','o','MarkerSize',4);
    hold on
    set(h,'LineStyle','none','MarkerEdgeColor','none')
    subplot(sub_h(iS));iS = iS+1;
    h = plot(wf_prop1(cell_inx),wf_prop4(cell_inx),'MarkerFaceColor',mk_color(co,:),'Marker','o','MarkerSize',4);
    hold on
    set(h,'LineStyle','none','MarkerEdgeColor','none')
    subplot(sub_h(iS));iS = iS+1;
    h = plot(wf_prop2(cell_inx),wf_prop3(cell_inx),'MarkerFaceColor',mk_color(co,:),'Marker','o','MarkerSize',4);
    hold on
    set(h,'LineStyle','none','MarkerEdgeColor','none')
    subplot(sub_h(iS));iS = iS+1;
    h = plot(wf_prop2(cell_inx),wf_prop4(cell_inx),'MarkerFaceColor',mk_color(co,:),'Marker','o','MarkerSize',4);
    hold on
    set(h,'LineStyle','none','MarkerEdgeColor','none')
    subplot(sub_h(iS));iS = iS+1;
    h = plot(wf_prop3(cell_inx),wf_prop4(cell_inx),'MarkerFaceColor',mk_color(co,:),'Marker','o','MarkerSize',4);
    hold on
    set(h,'LineStyle','none','MarkerEdgeColor','none')
end

% Plot light evoked spikes.
lspk = [];
spk = loadcb(cellid,'Spikes');
for iT = 1:length(pon),
    t_lspk = spk(spk > pon(iT) & spk < pon(iT) + 0.002);
    %    t_lspk = all_spikes(find(all_spikes > pon(iT) & all_spikes < pon(iT) + 0.002));
    % take the first light evoked spike.
    %     if length(t_lspk) > 1,
    %         disp(1)
    %         t_lspk = t_lspk(1);
    %     end
    lspk  = [lspk ;t_lspk];
end
[junk,junk2,cell_inx] = intersect(lspk,all_spikes);

% overlay light evoked spikes on plot.
axes(ax1)
h = plot(wf_prop2(cell_inx),wf_prop4(cell_inx),'MarkerFaceColor',[0 0.45 1],'Marker','o','MarkerSize',6);
set(h,'LineStyle','none','MarkerEdgeColor','none')
l(1) = xlabel('Energy 2');
l(2) = ylabel('Energy 4');
setmyplot(gca,l)
box off
axis square
axis([500 9000 500 6000])
set(gca,'XTick',xlim,'YTick',ylim)

% Plot light evoked spikes on all plots.
figure(f2h)
iS = 1;
subplot(sub_h(iS));iS = iS+1;
h = plot(wf_prop1(cell_inx),wf_prop2(cell_inx),'MarkerFaceColor',[0 0.45 1],'Marker','o','MarkerSize',6);
set(h,'LineStyle','none','MarkerEdgeColor','none')
subplot(sub_h(iS));iS = iS+1;
h = plot(wf_prop1(cell_inx),wf_prop3(cell_inx),'MarkerFaceColor',[0 0.45 1],'Marker','o','MarkerSize',6);
set(h,'LineStyle','none','MarkerEdgeColor','none')
subplot(sub_h(iS));iS = iS+1;
h = plot(wf_prop1(cell_inx),wf_prop4(cell_inx),'MarkerFaceColor',[0 0.45 1],'Marker','o','MarkerSize',6);
set(h,'LineStyle','none','MarkerEdgeColor','none')
subplot(sub_h(iS));iS = iS+1;
h = plot(wf_prop2(cell_inx),wf_prop3(cell_inx),'MarkerFaceColor',[0 0.45 1],'Marker','o','MarkerSize',6);
set(h,'LineStyle','none','MarkerEdgeColor','none')
subplot(sub_h(iS));iS = iS+1;
h = plot(wf_prop2(cell_inx),wf_prop4(cell_inx),'MarkerFaceColor',[0 0.45 1],'Marker','o','MarkerSize',6);
set(h,'LineStyle','none','MarkerEdgeColor','none')
subplot(sub_h(iS));iS = iS+1;
h = plot(wf_prop3(cell_inx),wf_prop4(cell_inx),'MarkerFaceColor',[0 0.45 1],'Marker','o','MarkerSize',6);
set(h,'LineStyle','none','MarkerEdgeColor','none')
