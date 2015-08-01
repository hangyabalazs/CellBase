function valid_channels = check_channel_validity(cellid)
%CHECK_CHANNEL_VALIDITY   Check tetrode channels.
%   VALID_CHANNELS = CHECK_CHANNEL_VALIDITY(CELLID) loads 'Energy' property
%   for all spikes on tetrode of CELLID and returns a 1-by-4 logical array
%   that is true for the tetrode channels that have a mean Energy higher
%   than an empirically decided threshold value.
%
%   See also LRATIO.

%   Edit log: BH 5/9/12

% Load feature data for tetrode.
[r,s,t,u] = cellid2tags(cellid);

% Load Energy
propfn = [getpref('cellbase','cell_pattern') num2str(t) '_Energy'];
propfn_path = [cellid2fnames(cellid,'sess') filesep 'FD'];
if ~isdir(propfn_path)
    propfn_path = cellid2fnames(cellid,'sess');
end
propfn_full = [propfn_path filesep propfn];
wf_prop = load([propfn_full '.fd'],'-mat');
wf_prop = wf_prop.FeatureData;

% Channel is valid if Energy is not zero
mne = mean(wf_prop,1);
valid_channels = abs(mne) > 150;