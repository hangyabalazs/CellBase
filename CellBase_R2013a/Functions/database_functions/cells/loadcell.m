function stimes = loadcell(cellid,varargin)
%LOADCELL   Load spike & event times of a cell.
%   TS = LOADCELL(CELLID) loads spike and trial event times into the base 
%   workspace and returns spike times in TS.
%
%   STIMES = LOADCELL(CELLID,'CALLER') loads the variables into the caller
%   functions workspace.
%
%   To make analysis routines generally usable, LOADCELL allows passing
%   variables other than a cell ID. If a numeric vector is passed instead of
%   a cell ID, then it is simply returned as an output without error. This
%   allows analysis routines to be called for a spike time vector without
%   using CellBase.
%
%   See also LOADCB.

%   Edit log: BH 7/6/12

% Load CellBase
load(getpref('cellbase','fname'));

% Default workspace to load variables into
WS = 'base';

% Use user-specified workspace if passed
if nargin > 1
    if ismember(lower(varargin{1}),{'base','caller'})
        WS = lower(varargin{1});
    else
        disp('LOADCELL: wrong argument');
        return
    end
end

% Load spike and event times
if ~isnumeric(cellid)
    [fname_spikes, fname_events] = cellid2fnames(cellid);
    if ~exist(fname_spikes,'file')   % file not found
        disp(sprintf('File %s not found.',fname_spikes))
        stimes = 0;
        return
    end
    load_spikes = sprintf('TS=load(''%s'');',fname_spikes);
    load_events = sprintf('load(''%s'');',fname_events);
    evalin(WS,load_events);
    eval(load_spikes);
    stimes = [TS.TS];
else    % isnumeric(cellid)
    stimes = cellid;   % pass on the input
end