function  [number,nontetrodepartners] = nontetrodepairs(cellid)
%NONTETRODEPAIRS   Find cells from different tetrodes in a session.
%   [NUMBER,NONTETRODEPARTNERS] = NONTETRODEPAIRS(CELLID) returns the total 
%   number of units on other tetrodes as well as the list of
%   cell IDs from other tetrodes. 
%
%   See also CELLID2TAGS and TETRODEPAIRS.

%   Edit log: BH 1/24/2012

% Input argument check
error(nargchk(1,1,nargin))

% Load CellBase
load(getpref('cellbase','fname'));
allcells = CELLIDLIST;

% Find tetrode pairs
[ratname session tetrode unit] = cellid2tags(cellid);
searchfor = sprintf('%s_%s_%d',ratname,session,tetrode);
mypairs_tp = strmatch(searchfor,allcells);
searchfor = sprintf('%s_%s_%d',ratname,session);
mypairs_sp = strmatch(searchfor,allcells);
nontetrodepartners = allcells(setdiff(mypairs_sp,mypairs_tp));
number = length(nontetrodepartners);