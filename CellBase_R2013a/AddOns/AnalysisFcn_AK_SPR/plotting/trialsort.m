function [ind, part1, part2] = trialsort(s,partition1,partition2,trigger_event,sort_event,vt);
%
% ugly as hell
%
%partition1 = 'stimulus'
%partition2 = 'outcome'
%trigger_event = 'WaterPokeIn';
%sort_event = 'WaterPokeOut';

if iscellstr(sort_event)
     sort_event = sort_event{1};
end

i = 1;
switch lower(partition1)
    case 'outcome'
       TOSORT(i,1:length(s.Correct(vt))) = NaN;     %include all
       TOSORT(i,find(s.Correct(vt))) = 1;
       TOSORT(i,find(s.Error(vt))) = 2;
   case 'direction'
       TOSORT(i,1:length(s.Left(vt))) = NaN;     %include all
       TOSORT(i,find(s.Left(vt))) = 1;
       TOSORT(i,find(s.Right(vt))) = 2;
   case 'stimulus'
       TOSORT(i,:) = s.OdorValveID(vt);
   otherwise
       i = i-1;
end
i=i+1;
switch lower(partition2)
    case 'outcome'
       TOSORT(i,1:length(s.Correct(vt))) = NaN;     %include all
       TOSORT(i,find(s.Correct(vt))) = 1;
       TOSORT(i,find(s.Error(vt))) = 2;
   case 'direction'
       TOSORT(i,1:length(s.Left(vt))) = NaN;     %include all
       TOSORT(i,find(s.Left(vt))) = 1;
       TOSORT(i,find(s.Right(vt))) = 2;
   case 'stimulus'
       TOSORT(i,:) = s.OdorValveID(vt);
   otherwise
       i = i-1;
end
i=i+1;

if ~isempty(sort_event)
    durations =  eval(['s.' sort_event ' - s.' trigger_event]);
    TOSORT(i,:) = durations(vt);
else
   TOSORT(i,:) = 1:length(s.Correct(vt));
end

[x, ind] = sortrows(TOSORT');

part1=[];
part2=[];

if size(TOSORT,1) > 1
  x(end+1,:) = -2;
  junk = x;
  junk(isnan(junk))=-1;
 
  edge=find(diff(junk(:,1)));
  part1 = [[1;edge(1:end-1)] edge];
  for i=1:size(part1,1);
    part1(i,3) = mean(x(part1(i,1):part1(i,2),1));
  end
  part1(:,1:2)=abs(part1(:,1:2)-(length(x)+1));
  
  if size(TOSORT,1) > 2
    
    edge=find(diff(junk(:,2)));
    part2 = [[1;edge(1:end-1)] edge];  
    for i=1:size(part2,1);
     part2(i,3) = mean(x(part2(i,1):part2(i,2),2));  
    end
    
   part2(:,1:2)=abs(part2(:,1:2)-(length(x)+1));
  end % 2
end % 1

 