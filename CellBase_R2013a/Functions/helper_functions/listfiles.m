function  list = listfiles(varargin)
%LISTFILES   List files in a directory.
%   L = LISTFILES lists all files in the current directory.
%
%   L = LISTFILES(DR) lists all files in the current DR.
%
%   L = LISTFILES(DR,STR) lists all files in the current DR that contains
%   the string STR.
%
%   See also LISTDIR.

%   Edit log: BH 3/21/11

% Input argument check
if nargin == 0
    mypath = pwd;
else
    mypath = varargin{1};
end

% List directory
l = dir(mypath);
list = {l(~[l.isdir]).name};

% Check string input in file names and restrict output list
if nargin == 2
    pat = char(varargin{2});
    s = regexp(list,pat);
    matched = ~cellfun(@isempty,s);
    list = list(matched);
end