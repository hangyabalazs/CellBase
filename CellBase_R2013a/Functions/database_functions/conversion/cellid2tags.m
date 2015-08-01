function   [ratname,session,tetrode,unit] = cellid2tags(cellid)
%CELLID2TAGS   Convert cell IDs to separate tags.
%   [RATNAME,SESSION,TETRODE,UNIT] = CELLID2TAGS(CELLID) parses CELLID
%   string to RATNAME, SESSION, TETRODE and UNIT strings.
%
%   See also CELLID2FNAMES.

%   Edit log: BH 6/23/11

% Convert cell ID
if iscellid(cellid)   % cell ID
    cellid = char(cellid);
    [ratname,remain] = strtok(cellid,'_');
    [session,remain] = strtok(remain(2:end),'_');
    [tetrodeunit] = sscanf(remain(2:end),'%d.%d');
    tetrode = tetrodeunit(1);
    unit = tetrodeunit(2);
elseif issessionid(cellid)   % session ID
    ratname = cellid{1};
    session = cellid{2};
    tetrode = [];
    unit = [];
else
    error('Unknown cell ID format.')
end