function BatchSessionClust(sessionpath)

% Make a PDF which contains cluster information of all units
% Saved in the session dir

if nargin==0;
[sessionpath]=uigetdir ('specify a session path');
end

cd(sessionpath);
LIST=dir([sessionpath, '\' 'TT*.clusters']);
ClustNames=char(LIST.name);

for it=1:size(ClustNames, 1)
% [clname]=strtok(ClustNames(it, :), '.clusters');
clname=ClustNames(it, :);
ClusterChecker(clname,sessionpath)
end