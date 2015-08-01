function nlxpos2mat_v4(varargin)
% This function converts the position information from nvt Neuralynx Video Tracker format ...
% into Mat files
% this replaces the MatlabRead function. It can convert a specified
% continuous channel from Neuralynx to matlab. 
% INPUTS:
% path of the data directory
% Continuous channel no. i.e. 1, 2,etc or 'Events' to convert events file
% if varargin(2) = 'Events' then read Events.nev into mat file
% If nothing is specified it converts all files from chan 4 to 8.
% SPR 2008-10-30


if nargin < 1,  
    [fn,path] = uigetfile('C:\Users\Sachin\Documents\','Select a file to transform');
    cd(path(1:end-1));
    pname = [path '\' fn];
else
    pname = varargin{1};
end

param2 = [1 1 1 1 1 1]; % Field Selection
param3 = 1; % Extract Header
param4 = 1; % Extract Mode
param5 = []; % Extraction Mode Array

[TimeStamps, ExtractedX, ExtractedY, ExtractedAngle, Targets, Points, NlxHeader] = Nlx2MatVT_v4( pname,param2,param3,param4,param5);
% [TimeStamps, ExtractedX, ExtractedY, ExtractedAngle, Targets, Points] = Nlx2MatVT_v4( Filename, FieldSelection, ExtractHeader, ExtractMode, ModeArray );

save([path '\' fn(1:end-4)],'TS','ScNumbers', 'CellNumbers', 'Params', 'DataPoints','NlxHeader');