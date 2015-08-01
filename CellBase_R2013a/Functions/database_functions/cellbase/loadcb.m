function  varargout = loadcb(varargin)
%LOADCB   Load the CellBase files.
%   LOADCB loads CellBase (CELLIDLIST, ANALYSES, TheMatrix) into the caller
%   function's workspace.
%
%   OUT = LOADCB(CELLID,FILETYPE) loads the specified file (FILETYPE) of
%   a given cell (CELLID) to OUT, where filetype may be 'Spikes', 
%   'TrialEvents', 'StimEvents', 'STIMSPIKES', 'EVENTSPIKES', 'Continuous'
%   or 'Waveforms' (see also CELLID2FNAMES).
%
%   See also INITCB and CELLID2FNAMES.

% Edit log: AK 3/04, 11/06, 4/10; BH 7/6/12

% Without input arguments, load the entire database
if nargin == 0
    if ~isempty(getpref('cellbase')) 
        evalin('caller','load(getpref(''cellbase'',''fname''))');
    else
        disp('LOADCB: CellBase does not exist.');
    end
else 
    
    % Load variables for a specified cell
    cellid = varargin{1};   % we must have a cellid
    
    % Filetype to load
    if nargin > 1
        filetype = varargin{2};
    else
        filetype = 'Spikes';   % the default is to load spikes
    end
    
    % Filename to load
    fname = cellid2fnames(cellid,filetype);
    TIMEFACTOR = getpref('cellbase','timefactor');    % scaling factor to convert spike times into seconds
    
    if nargout == 1     % if we are loading into a variable
        if strncmpi(filetype,'Spikes',5)
            x = load(fname);
            varargout{1} = x.TS*TIMEFACTOR;  % if we are loading spikes, then multiply with conversion factor
        elseif strncmpi(filetype,'Waveforms',4)
            SpikeTimes = loadcb(cellid,'Spikes');  % load stimulus spikes (prealigned)
            
            % Load waveform data (Ntt file)
            Nttfile = cellid2fnames(cellid,'ntt');
            TIMEFACTOR = getpref('cellbase','timefactor');    % scaling factor to convert spike times into seconds
            [all_spikes all_waves] = LoadTT_NeuralynxNT(Nttfile);
            [junk junk2 evoked_inx] = intersect(SpikeTimes,all_spikes*TIMEFACTOR);
            if ~isequal(junk,SpikeTimes)   % internal check for spike times
                error('loadcb:SpikeTimeMismatch','Mismatch between extracted spike times and Ntt time stamps.')
            end
            varargout{1} = all_waves(evoked_inx,:,:);   % selected waveforms
        else
            varargout{1} = load(fname);
        end
    else     % otherwise load into the workspace
        if strncmpi(filetype,'Spikes',5) || strncmpi(filetype,'Waveforms',4)
            error('loadcb:noOutputArg','Specify output argument for LOADCB.')
        else
            evalin('caller',sprintf('load(''%s'')',fname));
        end
    end
end