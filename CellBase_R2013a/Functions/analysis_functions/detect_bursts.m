function Bstimes = detect_bursts(stimes)
%DETECT_BURSTS   Burst detection.
%   BSTIMES = DETECT_BURSTS(STIMES) finds bursts in spike time array
%   STIMES (in seconds). Bursts are defined as spikes with ISIs
%   corresponding to at least 120 Hz.
%
%   See also VIEWCELL2B.

%   Edit log: BH 6/23/11

% Detect bursts
lens = length(stimes);
Bstimes = cell(1,lens);
ISIstimes = cell(1,lens);
for i = 1:lens
    for k = 1:length(stimes{1,i}) - 1
        ISIstimes{1,i}(k) = 1 / (stimes{1,i}(k+1) - stimes{1,i}(k));
        if ISIstimes{1,i}(k) >= 120;
            Bstimes{1,i}(k) = stimes{1,i}(k);
        else
            Bstimes{1,i}(k) = NaN;
        end
    end
end