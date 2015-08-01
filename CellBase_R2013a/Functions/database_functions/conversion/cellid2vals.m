 function [ratnum ses tetrode unit] = cellid2vals(cellid)
%CELLID2VALS   Convert cell IDs to numbers.
%   [R S T U] = CELLID2VALS(CELLID) converts cell ID to 4 numbers: serial
%   position for the rat (R), date number in yymmdd format (S), tetrode (T)
%   and unit (U) number.
%
%   See also CELLID2TAGS.

%   Edit log: BH 6/25/12, 7/5/12

% Get tags
[ratname,session,tetrode,unit] = cellid2tags(cellid);

% Convert ratname to a serial number
rats = sort(listdir(getpref('cellbase','datapath')));
ratnum = strmatch(ratname,rats);

% Keep only the date part of session name
ses = str2double(session(1:end-1));