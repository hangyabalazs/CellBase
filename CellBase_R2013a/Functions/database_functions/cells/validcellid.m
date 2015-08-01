function yn = validcellid(cellid,varargin)
%VALIDCELLID   Check if CellBase cell ID is valid.
%   YN = VALIDCELLID(CELLID) returns 1 if the referred file exists, 0
%   if it doesn't and -1 if the cell ID couldn't be parsed.
%
%   YN = VALIDCELLID(CELLID,'LIST') performs the search in CellBase 
%   CELLIDLIST.
%
%   See also CELLID2FNAMES and ISCELLID.

%   Edit log: BH 6/23/11, 5/4/12

% Search
if nargin > 1 &&  strcmpi(varargin{1},'list')   % search in CELLIDLIST     
    load(getpref('cellbase','fname'),'CELLIDLIST');
    yn = ~isempty(strmatch(char(cellid),CELLIDLIST));
else
    try
        fname_spikes = cellid2fnames(cellid);   % without 'list' option
    catch
         yn = -1; % error
         return
    end
    yn = min(exist(fname_spikes,'file'),1);
end