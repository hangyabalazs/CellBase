function  list = listtag(xstr)
%LISTTAG   List CellBase tags.
%   L = LISTTAG(STR) returns a cell array containing the specified tag from
%   CellBase. Options for STR:
%   'properties', 'analysis', 'cell', 'rat'/'animal', 'session', 'tetrode'.
%   'properties' and 'analysis' are searched in ANALYSES, while others are
%   searched in CELLIDLIST. The 'allsession' option includes sessions
%   without cells (e.g. behavior only) by searching the directory tree.
%
%   See also LISTDIR and LISTFILES.

%   Edit log: BH 3/21/11, 5/3/12, 5/20/14

% Load cellbase
global CELLIDLIST ANALYSES TheMatrix
if isempty(CELLIDLIST)
    load(getpref('cellbase','fname'));
end

% Get tagged list
xstr = lower(char(xstr));
list = {};
if strncmp(xstr,'prop',4)
    n = 1;
    for i = 1:length(ANALYSES)
        for j = 1:length(ANALYSES(i).propnames)
            list{n} = char(ANALYSES(i).propnames(j));
            n = n + 1;
        end
    end
elseif strncmp(xstr,'anal',4)
    for i = 1:length(ANALYSES)
        list{i} = func2str(ANALYSES(i).funhandle);
    end
elseif strncmp(xstr,'rat',3) || strncmp(xstr,'animal',3) || strncmp(xstr,'mouse',3) ...
        || strncmp(xstr,'mice',3)
    clist = char(CELLIDLIST);
    for i = 1:length(clist)
        [rat,remain] = strtok(clist(i,:),'_');
        list{i} = rat;
    end
    list = unique(list);
elseif strncmp(xstr,'ses',3)
    clist = char(CELLIDLIST);
    for i = 1:length(clist)
        [rat,remain] = strtok(clist(i,:),'_');
        [session,remain] = strtok(remain(2:end),'_');
        list{i,1} = rat(:)';
        list{i,2} = session(:)';
    end
    list = unique_cell(list);
elseif strncmp(xstr,'allses',6)
    clist = dir(getpref('cellbase','datapath'));
    clist = clist(3:end);
    idr = [clist.isdir];
    clist = {clist(idr).name};
    list = cell(0,2);
    for i = 1:length(clist)
        clist2 = dir(fullfile(getpref('cellbase','datapath'),clist{i}));
        clist2 = clist2(3:end);
        idr = [clist2.isdir];
        clist2 = {clist2(idr).name};
        nums = length(clist2);
        rat = repmat(clist(i),nums,1);
        session = clist2';
        ll = size(list,1);
        list(ll+1:ll+nums,1) = rat;
        list(ll+1:ll+nums,2) = session;
    end
    list = unique_cell(list);
elseif strncmp(xstr,'tetrode',3)
    clist = char(CELLIDLIST);
    for i = 1:length(clist)
        [rat,remain] = strtok(clist(i,:),'_');
        [session,remain] = strtok(remain(2:end),'_');
        [tetrode,remain] = strtok(remain(2:end),'.');
        list{i,1} = rat;
        list{i,2} = session;
        list{i,3} = tetrode;
    end
    list = unique_cell(list);
elseif strncmp(xstr,'cell',3)
    list = CELLIDLIST;
else
    disp('LISTTAG: unrecognized option.')
end