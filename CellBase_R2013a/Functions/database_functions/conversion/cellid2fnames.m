function [fname_spikes,fname_events] = cellid2fnames(cellid,filename,CSC_chan)
%CELLID2FNAMES   Convert cell IDs to data file names.
%   [FNAME_SPIKES, FNAME_EVENTS] = CELLID2FNAMES(CELLID) returns spike data
%   and event filenames for CELLID.
%
%   [FNAME_SPIKES, FNAME_EVENTS] = CELLID2FNAMES(CELLID,FILENAME) uses
%   FILENAME input argument to genearte FNAME_SPIKES. Specials cases for
%   FILENAME are handled differentially:
%       'TrialEvent', behavioral session name ('TrialEvents2.mat') is returned
%       'StimEvent', stimulation filename is returned
%       'Session', directory path for the session
%       'Position', 'POSITION' string is used
%       'Spikes', 'cell_pattern' preference is used to generate .mat filename
%       'tfile', 'cell_pattern' preference is used to generate .t filename
%       'Ntt', .Ntt filename is returned
%       'wv', 'cell_pattern' preference is used to return *wv.mat filename
%       'quality', *ClusterQual.mat filename is returned
%       'Waveform', waveform data filename is returned
%       'Events', 'EVENTS' string is used
%       'Eventspikes', 'EVENTSPIKES' string is used
%       'Stimspikes', 'STIMSPIKES' string is used
%       'Continuous', CSCx_x.mat or CSCx.mat filename is returned
%       'Evewaves', EVENTSPIKESxcx.mat filename is returned
%       'stimwaves', STIMWAVESxcx.mat filename is returned.
%   For the 'Continuous', 'Evewaves' and 'stimwaves' options, LFP channel
%   names should be provided as a third input argument, using
%   [FNAME_SPIKES, FNAME_EVENTS] = CELLID2FNAMES(CELLID,FILENAME,CSC_CHAN)
%   syntax.
%
%   See also CELLID2TAGS.
    
%   Edit log: ZFM 10/7/04, AK 11/06, AK 4/10, SPR 07/2010, BH 6/23/11

% Get CellBase preferences
cellbase_datapath = getpref('cellbase','datapath');
session_fname = getpref('cellbase','session_filename');
if ispref('cellbase','StimEvents_filename')
	stim_fname = getpref('cellbase','StimEvents_filename');
end
cellbase_cell_pattern = getpref('cellbase','cell_pattern');
continuous_channel = 'CSC';

% Get tags
[ratname,session,tetrode,unit] = cellid2tags(cellid);

% Create unit filename
tetrodeunit = sprintf('%s%d_%d.mat',cellbase_cell_pattern,tetrode,unit);

% Create names
if nargin < 2   % if filename was not specified
    fname_spikes = fullfile(cellbase_datapath,ratname,session,tetrodeunit);
    fname_events = fullfile(cellbase_datapath,ratname,session,session_fname);
else
    % not really spikes, but whatever you specified
    % create unit filename
    if strncmpi(filename,'TrialEvent',10)
        fname_unit = session_fname;     %'TrialEvents2.mat';
    elseif strncmpi(filename,'StimEvent',9)
        fname_unit = stim_fname;
    elseif strncmpi(filename,'Session',3)
        fname_unit = '';
    elseif strncmpi(filename,'Position',3),
        fname_unit='POSITION';
    elseif strncmpi(filename,'Spikes',5)
        fname_unit = sprintf('%s%d_%d.mat',cellbase_cell_pattern,tetrode,unit);
    elseif strncmpi(filename,'tfile',5)
        fname_unit = sprintf('%s%d_%d.t',cellbase_cell_pattern,tetrode,unit);
    elseif strncmpi(filename,'Ntt',3)
        fname_unit = sprintf('TT%d.ntt',tetrode);
    elseif strncmpi(filename,'wv',2)
        fname_unit = sprintf('%s%d_%d-wv.mat',cellbase_cell_pattern,tetrode,unit);
    elseif strncmpi(filename,'quality',4)
        fname_unit = sprintf('%s%d_%d-ClusterQual.mat',cellbase_cell_pattern,tetrode,unit);
    elseif strncmpi(filename,'Waveform',8)
        fname_unit = sprintf('WAVEFORMDATA%d_%d.mat',tetrode,unit);
    elseif strmatch(filename,'Events','exact')
        fname_unit = 'EVENTS';   % SPR 09/12/24
    elseif strncmpi(filename,'Eventspikes',8)
        fname_unit = sprintf('EVENTSPIKES%d_%d.mat',tetrode,unit);
    elseif strncmpi(filename,'Stimspikes',6)
        fname_unit = sprintf('STIMSPIKES%d_%d.mat',tetrode,unit);        
    elseif strncmpi(filename,'Continuous',4) || strncmpi(filename,'CSC',3)
        if nargin < 3   % you have to specify which LFP you want
            warning('cellid2fnames:inputargMissing','Specify LFP channel with a [1 X 2] vector using CSC filename convention')
            tetrode = 0;
            chan = 1;
        else
            tetrode = CSC_chan(1);
            if length(CSC_chan) > 1
                chan = CSC_chan(2);
            end
        end
        if isequal(length(CSC_chan),1)
            fname_unit = sprintf('%s%d.mat',continuous_channel,tetrode); %CSCx.mat
        else
            fname_unit = sprintf('%s%dc%d.mat',continuous_channel,tetrode,chan); %CSCx_x.mat
        end
    elseif strncmpi(filename,'Evewaves',4),
        if nargin < 3   % you have to specify which LFP you want
            warning('cellid2fnames:inputargMissing','Specify LFP channel with a [1 X 2] vector using CSC filename convention')
            tetrode = 0;
            chan = 1;
        else
            tetrode = CSC_chan(1);
            chan = CSC_chan(2);
        end
        fname_unit = sprintf('%s%dc%d.mat',filename,tetrode,chan);    % EVENTSPIKESxcx.mat
    elseif strncmpi(filename,'stimwaves',6),
        if nargin < 3   % you have to specify which LFP you want
            warning('cellid2fnames:inputargMissing','Specify LFP channel with a [1 X 2] vector using CSC filename convention')
            tetrode = 0;
            chan = 1;
        else
            tetrode = CSC_chan(1);
            chan = CSC_chan(2);
        end
        fname_unit = sprintf('%s%dc%d.mat',filename,tetrode,chan);    % STIMWAVESxcx.mat
    else
        fname_unit = sprintf('%s%d_%d.mat',filename,tetrode,unit);
    end
    fname_spikes = fullfile(cellbase_datapath,ratname,session,fname_unit);
end