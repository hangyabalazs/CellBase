%%

mousedir = 'f:\HDB_cellbase\n079\';

mousename = mousedir(end-4:end-1);
mousename2 = ['nb0' mousedir(end-2:end-1)];
subdirs = dir(mousedir);
subdirs = {subdirs.name};
subdirs = subdirs(3:end);

DirNum = length(subdirs);
CHK = false(1,DirNum);
for k = 1:DirNum
    cdir = subdirs{k};
    inpdir = fullfile(getpref('cellbase','datapath'),mousename,cdir);
    fout = checknlxconfig(inpdir);
    if isequal(fout,-1)
        disp([inpdir ': Config file not found.'])
        CHK(k) = true;
    elseif isequal(fout,-2)
        disp([inpdir ': Previous config file was used.'])
        CHK(k) = true;
    else
        CHK(k) = strcmp(mousename2,fout(1:5));
        if ~CHK(k)
            disp([inpdir ': ' fout ' config file was used.'])
        end
    end
end

disp(CHK)
if ~all(CHK)
    error('Screwed up.')
end