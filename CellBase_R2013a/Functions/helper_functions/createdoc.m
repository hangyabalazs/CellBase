function createdoc(fn,resdir)
%CREATEDOC   Create documentation html file for functions.
%   CREATEDOC(FN,DR) writes an html documentation file for the function FN
%   to the directory DR. It parses the help of FN, which should be
%   formulated according to Matlab conventions. Custom Matlab R2010a html
%   design is used.
%
%   See also HELP and DOC.

%   Balazs Hangya
%   Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   New York 11724, USA
%   balazs.cshl@gmail.com

% Parse filename
[pathname,filename,extension] = fileparts(which(fn));

% Parse help
[nm dsc syntx hlpbody hlpbody_cell seealso] = parsehelp(fn);
for hh = 1:length(hlpbody_cell)
    thlp = hlpbody_cell{hh};
    try
        wrds2 = strread(thlp,'%s','delimiter',' ');   % add tt commands
        sps = [regexp(thlp,' ') length(thlp)];
        pap1 = [false; strcmp(wrds2,cellfun(@upper,wrds2,'UniformOutput',false)); false];
        pap2 = [false; ~isbetweentts(thlp); false];
        pap3 = [false; isnan(str2double((wrds2))); false];
        allcap = pap1 & pap2 & pap3;
        da = diff(allcap);
        fda = find(da==1);
        fda2 = find(da==-1);
        cntr = 1;
        while ~isempty(fda) && cntr < 30
            if fda(1) == 1
                inx = 0;
            else
                inx = sps(fda(1)-1);
            end
            thlp = [thlp(1:inx) '<tt>' ...
                lower(thlp(inx+1:sps(fda2(1)-1)-1))...
                '</tt>' thlp(sps(fda2(1)-1):end)];
            sps = [regexp(thlp,' ') length(thlp)];
            wrds2 = strread(thlp,'%s','delimiter',' ');
            pap1 = [false; strcmp(wrds2,cellfun(@upper,wrds2,'UniformOutput',false)); false];
            pap2 = [false; ~isbetweentts(thlp); false];
            pap3 = [false; isnan(str2double((wrds2))); false];
            allcap = pap1 & pap2 & pap3;
            da = diff(allcap);
            fda = find(da==1);
            fda2 = find(da==-1);
            cntr = cntr + 1;   % counter to prevent infinite loop
        end
    catch ME
        disp('Adding tt commands failed.')
        disp(ME.message)
    end
    hlpbody_cell{hh} = thlp;
end

% Write html file
try
    fp = filesep;
    fid = fopen([resdir fp filename '.html'],'w');
    writeheader(fid,nm)    % write html header
    fprintf(fid,'%s\n','<body bgcolor="#FFFFFF" text="#000001">');
    fprintf(fid,'%s\n',['<h1 align="left">' nm '</h1>']);
    fprintf(fid,'%s\n',['<align="left">' dsc]);
    if ~isempty(syntx)
        fprintf(fid,'%s\n','<h2>Syntax</h2>');
        fprintf(fid,'%s','<p><blockquote>');
        for k = 1:length(syntx)
            fprintf(fid,'%s\n',['<tt>' lower(syntx{k}) '</tt><br>']);
        end
        fprintf(fid,'%s\n','</blockquote></p>');
    end
    fprintf(fid,'%s\n','<h2>Description</h2>');
    for k = 1:length(hlpbody_cell)
        thlp = hlpbody_cell{k};
        fprintf(fid,'%s\n',['<p>' thlp '</p>']);
    end
    fprintf(fid,'%s\n','<h2>See Also</h2>');
    lens = length(seealso);
    fprintf(fid,'%s\n','<p>');
    for k = 1:lens
        if isempty(findstr(which(seealso{k}),'toolbox'))
            fprintf(fid,'%s',['<a href="' seealso{k} '.html"><tt>' seealso{k} '</tt></a>']);
        else
            fprintf(fid,'%s',['<a href="matlab:doc ' seealso{k} '"><tt>' seealso{k} '</tt></a>']);
        end
        if k < lens
            fprintf(fid,'%s',', ');
        else
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'%s','</p>');
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'%s\n','</body>');
    fprintf(fid,'%s\n','</html>');
    fclose(fid);
catch ME
    fclose(fid);
    rethrow(ME)
end

% -------------------------------------------------------------------------
function [nm dsc syntx hlpbody hlpbody_cell seealso] = parsehelp(fn)

% Parse help
hlp = help(fn);
rh = regexp(hlp,'\n');
pdsc = hlp(1:rh(1)-1);
fn = find(pdsc~=' ',1,'first');  % skip spaces
pdsc = pdsc(fn:end);
rh2 = regexp(pdsc,' ');
nm = lower(pdsc(1:rh2(1)-1));   % function name
pdsc2 = pdsc(rh2(1)+1:end);
fn = find(pdsc2~=' ',1,'first');
dsc = pdsc2(fn:end);  % function description
phb = hlp(rh(1)+1:rh(end-1)-1);
da = diff([0 phb==' ' 0]);   % eliminate multiple places
fda = find(da==1);
fda2 = find(da==-1);
spsl = max(fda2-fda);
for k = spsl:-1:1
    sps = repmat(' ',1,k);
    phb = regexprep(phb,sps,' ');
end
fn = find(phb~=' ',1,'first');  % skip initial spaces
hlpbody = phb(fn:end);
ffn = [1 strfind(hlpbody,char([10 32 10])) length(hlpbody)];   % make cell array of paragraphs
lenffn = length(ffn) - 1;
hlpbody_cell = cell(1,lenffn);
for k = 1:lenffn
    hlpbody_cell{k} = strtrim(hlpbody(ffn(k):ffn(k+1)));
end
wrds = strread(hlpbody,'%s','delimiter',' ');
eqs = find(strcmp(wrds,'='));
leneqs = length(eqs);
syntx = cell(1,leneqs);
for k = 1:leneqs   % syntax
    str = [lower(wrds{eqs(k)-1}) ' = ' lower(wrds{eqs(k)+1})];
    syntx{k} = str;
end
psa = hlp(rh(end-1)+1:end);
wrds2 = strread(psa,'%s','delimiter',' .,');
allcap = find(strcmp(wrds2,cellfun(@upper,wrds2,'UniformOutput',false)));
seealso = cellfun(@lower,wrds2(allcap),'UniformOutput',false);   % see also

% -------------------------------------------------------------------------
function out = isbetweentts(str)

% Determine whether already in a tt-block
wrds = strread(str,'%s','delimiter',' ');
lenw = length(wrds);
out = false(lenw,1);
ttstart = find(~cellfun(@isempty,strfind(wrds,'<tt>')));
ttend = find(~cellfun(@isempty,strfind(wrds,'</tt>')));
for k = 1:lenw
    tts = ttstart(find(ttstart<k,1,'last'));
    if ~isempty(tts)
        tte = ttend(find(ttend>tts,1,'first'));
        if tte > k
            out(k) = true;
        end
    end
end

% -------------------------------------------------------------------------
function sps = regexp2(str,dlm)

% Match to either of the delimiters
sp = [];
for k = 1:length(dlm)
    sp = [sp findstr(str,dlm(k))];
end
sps = sort(unique(sp),'ascend');

% -------------------------------------------------------------------------
function writeheader(fid,nm)

% Write html header
fprintf(fid,'%s\n','<html>');
fprintf(fid,'%s\n','<head>');
fprintf(fid,'%s\n',['<title>' nm '</title>']);
fprintf(fid,'%s\n','<meta http-equiv="Content-Type" content="text/html; charset=utf-8">');
fprintf(fid,'%s\n','<style>');
fprintf(fid,'%s\n','body {');
fprintf(fid,'%s\n','  background-color: white;');
fprintf(fid,'%s\n','  margin:1px;');
fprintf(fid,'%s\n','}');
fprintf(fid,'%s\n','h1 {');
fprintf(fid,'%s\n','  color: #990000; ');
fprintf(fid,'%s\n','  font-size: 30;');
fprintf(fid,'%s\n','  font-weight: normal;');
fprintf(fid,'%s\n','  margin-top: 12px;');
fprintf(fid,'%s\n','  margin-bottom: 0px;');
fprintf(fid,'%s\n','}');
fprintf(fid,'%s\n','h2 {');
fprintf(fid,'%s\n','  color: #990000;');
fprintf(fid,'%s\n','  font-size: 20;');
fprintf(fid,'%s\n','  margin-top: 24px;');
fprintf(fid,'%s\n','  margin-bottom: 0px;');
fprintf(fid,'%s\n','}');
fprintf(fid,'%s\n','pre.codeinput {');
fprintf(fid,'%s\n','  margin-left: 30px;');
fprintf(fid,'%s\n','}');
fprintf(fid,'%s\n','bld {');
fprintf(fid,'%s\n','  font-weight: bold;');
fprintf(fid,'%s\n','}');
fprintf(fid,'%s\n','</style>');
fprintf(fid,'%s\n','</head>');
fprintf(fid,'\n');