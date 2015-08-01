session_cellids=unique_session_cells;
% newbasedir='/Users/ranades/Documents/work/Data/mPFC/mPFC cellbase lite';
% newbasedir='/Users/ranades/Documents/work/Data/mPFC/mPFC_cellbase_c_May2011';
% newbasedir='/Volumes/SACHDATA3/Data/mPFC/mPFC_cellbase';
newbasedir = 'c:\Balazs\_data\SOM_Sachin\mPFC_cellbase';
co=1;
for iSess=1:length(session_cellids),
    tic;
    cellid=session_cellids(iSess);
    parentdir=cellid2fnames(cellid,'Sess');
    [rat,sess]=cellid2tags(cellid);
    
    yo=strfind(parentdir,getpref('cellbase','datapath'));
    ratdir = [newbasedir filesep rat];
    if ~exist(ratdir,'dir')
        mkdir(ratdir)
    end
    sessdir = [ratdir filesep sess];
    if ~exist(sessdir,'dir')
        mkdir(sessdir)
    end
    destdir=[newbasedir parentdir(yo+length(getpref('cellbase','datapath')):end)];
    destdir=[sessdir];
    filelist=listfiles(parentdir);
    xclude_pattern='CSC';
    include_files=filelist(setdiff(1:length(filelist),strmatch(xclude_pattern,filelist)));
    for iF=1:length(include_files),
        if ~isempty(findstr('.mat',include_files{iF})),
            [success,msg,msgid]=copyfile([parentdir filesep include_files{iF}],destdir);
            if success == 0,
                failedfiles{co} = [parentdir filesep include_files{iF}];
                co = co+1;
                disp([parentdir filesep include_files{iF}]);
            end
        end
    end
    toc;
    disp([iSess length(session_cellids) cellid])
end