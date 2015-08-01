function  list = listdir(varargin)
%LISTDIR   List descending directories.
%   LIST = LISTDIR lists all the subdirectories in current directory,
%   except for '.' and '..'.
%
%   LIST = LISTDIR(DR) lists all the subdirectories in DR, except for '.'
%   and '..'.
%
%   See alse DIR and LISTFILES.

%   Edit log: BH 3/21/11

% Input argument check
if nargin == 0
    mypath = pwd;
else
    mypath = varargin{1};
end

% List subdirectories
d = dir(mypath);
list = {d([d.isdir]).name};
list = setdiff(list,{'.','..'});