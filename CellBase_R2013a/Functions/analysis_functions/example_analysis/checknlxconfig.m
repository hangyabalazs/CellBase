function fout = checknlxconfig(inpdir)

fnm = fullfile(inpdir,'CheetahLogFile.txt');
if ~exist(fnm,'file')
    fout = -1;
    return
end
contents = textread(fnm,'%s');
logfnd = strfind(contents,'Balazs_Configuration');
logfpos = find(cellfun(@(s)~isempty(s),logfnd));
if isempty(logfpos)
    logfnd = strfind(contents,'LastCheetah');
    logfpos = find(cellfun(@(s)~isempty(s),logfnd));
    if ~isempty(logfpos)   % previous setting was imported
        fout = -2;
        return
    end
end
logfile = contents{logfpos(1)};
[pth fout ext] = fileparts(logfile);

% keyboard