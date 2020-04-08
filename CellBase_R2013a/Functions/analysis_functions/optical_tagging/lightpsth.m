function lightpsth(sessionpath)
%LIGHTPSTH   Spike times aligned to photostimulation.
%   LIGHTPSTH(SESSIONPATH) calculates peri-event time histogram aligned to
%   the onset of photostimulation pulse trains for unclustered tetrode
%   data.
%
%   Saa also TAGGEDPROP.

%   Balazs Hangya
%   Institute of Experimental Medicine, Budapest, Hungary
%   hangya.ba;azs@koki.mta.hu

% List of TT files
dr = dir(sessionpath);
files = {dr.name};
TTpattern = getpref('cellbase','cell_pattern');   % find tetrode files with a cellbase-defined naming convention
TTfiles = strfind(files,TTpattern);
TTinx = regexp(files,[TTpattern '\d\.mat']);
TTinx = cellfun(@(s)~isempty(s),TTinx);
TTfiles = files(TTinx);

% Load photostimulation time stamps
fnm = fullfile(sessionpath,'lightevent.mat');
if ~exist(fnm,'file')
    convert_events(sessionpath);   % convert PulsePal events
end
load(fnm);

% Load spike times
NumTetrodes = length(TTfiles);
wn = [-20 100];   % PSTH window: -20 to 100 ms 
mwn = max(abs(wn));   % maximal lag
psth = nan(NumTetrodes,2*mwn+1);
legendstring = cell(1,NumTetrodes);
for iT = 1:NumTetrodes
    fnm = fullfile(sessionpath,TTfiles{iT});
    load(fnm)
    
    % Pseudo-trains
    mx = ceil(TimeStamps(end)*1000);
    pse = zeros(1,mx+1000);   % pseudo-event train, ms resolution
    pse(ceil(pulseon*1000)) = 1;
    psu = zeros(1,mx+1000);   % pseudo-spike train
    psu(ceil(TimeStamps*1000)) = 1;
    
    % PSTH
    [lpsth, lags] = xcorr(psu,pse,mwn);
    psth(iT,:) = lpsth;
    legendstring{iT} = ['TT' num2str(iT)];
end

% Plot
H = figure;
plot(lags,psth');
legend(TTfiles);
fnm = [sessionpath '\' 'light_psth.jpg'];   % save
saveas(H,fnm)
fnm = [sessionpath '\' 'light_psth.fig'];
saveas(H,fnm)
close(H)