function binraster = stimes2binraster(stimes,time,dt,ev_windows,window_margin)
%STIMES2BINRASTER   Calculate binraster from spike times.
%   BINRASTER = STIMES2BINRASTER(STIMES,TIME,DT,EV_WINDOWS,WINDOW_MARGIN)
%   calculates binrasters using spike times (STIMES) and time vector (TIME)
%   at a resolution of DT. Time window boundaries corresponding to the rows
%   of the binraster are specified by EV_WINDOWS. A margin of WINDOW_MARGIN
%   is added.
%
%   See also VIEWCELL2B and PLOT_RASTER2A.

%   Edit log: AK 7/1, BH 6/23/11, BH 3/11/12

% Input arguments check
NumTrials = length(stimes);
if ~exist('window_margin','var')
   window_margin = [0 0];
end
if exist('ev_windows','var')
    if size(ev_windows,1) ~= NumTrials
        ev_windows = ev_windows';
    end
    ev_windows = ev_windows + repmat(window_margin,NumTrials,1);
end
win_max = [1 length(time)];
    
% Spikes raster matrix
binraster = zeros(NumTrials,length(time));
for iTRIAL = 1:NumTrials
    all_spikes = stimes{iTRIAL};
    ok_spikes = all_spikes(all_spikes>time(1)&all_spikes<=time(end));
    ind_ok_spikes = round((ok_spikes-time(1))/dt) + 1;
    if ~isempty(ind_ok_spikes)
        binraster(iTRIAL,ind_ok_spikes) = 1;
    end
    if exist('ev_windows','var')
        win = ev_windows(iTRIAL,:);     
        ind_win = ceil((win-time(1))/dt);
        ind_win(isnan(win)) = win_max(isnan(win));   % replace NaNs with ends
        if ind_win(1) > ind_win(2)
            ind_win(2) = ind_win(1);
        end
        binraster(iTRIAL,[1:ind_win(1) max(ind_win(2),1):end]) = NaN;
    end
end   % iTRIAL