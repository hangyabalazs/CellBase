function [selts, seltsind, selisi] = extractSegSpikes(cellid,segs)
%EXTRACTSEGSPIKES   Extract spikes between specific timestamps.
%   [SELTS, SELTSIND, SELISI] = EXTRACTSEGSPIKES(CELLID,SEGS) extracts
%   timestamps for spikes of a cell (CELLID) that lie between start (first
%   row of SEGS) and end (second row of SEGS) points of given segments.
%   Output arguments:
%       SELTS: timestamps of the extracted spikes
%       SELTSIND: indices of the timestamps with respect to all spikes in
%           the cluster
%       SELISI: interspike intervals where more than one spike was found 
%           within a time segment; NaN is returned when there is only one
%           spike
%
%   EXTRACTSEGSPIKES(SPT,SEGS) overloads the function with numeric input in
%   the first argument, implemented as spike times.
%
%   See also FINDSEGS3 and ABS2RELTIMES.

%   Edit log: SPR 5/10, BH 4/25/12, 5/8/12

% Input argument check
prs = inputParser;
addRequired(prs,'cellid',@(s)iscellid(s)|isnumeric(s))
addRequired(prs,'segs',@isnumeric)

% Load spikes
if iscellid(cellid)
    spk = loadcb(cellid);
else
    spk = cellid;   % overload extractSegSpikes: implement first argument as spike times
end
spk = spk(:);

% Extract spikes for each segment
seltsind = [];
selisi = [];
for is = 1:size(segs,2)   % segment loop
    ladle = find(spk>segs(1,is)&spk<segs(2,is));
    seltsind = [seltsind; ladle]; %#ok<AGROW>
    if size(ladle,1)>1
        selisi = [selisi; diff(spk(ladle))]; %#ok<AGROW>
    else
        selisi = [selisi; NaN]; %#ok<AGROW>
    end
end
selts = spk(seltsind);