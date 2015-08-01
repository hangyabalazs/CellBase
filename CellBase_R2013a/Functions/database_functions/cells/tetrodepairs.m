function  [number,tetrodepartners] = tetrodepairs(cellid,includeself)
%TETRODEPAIRS   Find cells from the same tetrode in a session.
%   [NUMBER,TETRODEPARTNERS] = TETRODEPAIRS(CELLID) returns the total 
%   number of units on the same tetrode as CELLID as well as the list of
%   cell IDs from the same tetrode. 
%
%   TETRODEPAIRS(CELLID,INCLUDESELF) accepts a logical input argument
%   controling whether the unit for which the function was called is
%   included in the number and the list.
%
%   See also CELLID2TAGS.

%   Edit log: AK, SPR 5/10, BH 1/13/2012

% Input argument check
error(nargchk(1,2,nargin))
if nargin < 2
    includeself = true;     % default behavior: include the cellid for which the fcn was called
end

% Load CellBase
load(getpref('cellbase','fname'));
allcells = CELLIDLIST;

% Find tetrode pairs
[ratname session tetrode unit] = cellid2tags(cellid);
searchfor = sprintf('%s_%s_%d',ratname,session,tetrode);
mypairs = strmatch(searchfor,allcells);
tetrodepartners = allcells(mypairs);
number = length(mypairs);

% Drop 'self' if requested
if ~includeself
    number = number - 1;
    tetrodepartners = setdiff(tetrodepartners,cellid);
end