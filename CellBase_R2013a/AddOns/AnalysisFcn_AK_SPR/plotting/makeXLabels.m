function [xlabels, xvalues] = makeXLabels(fun_definitions,TAGS);


%AK 7/1

myXLabels = feval(fun_definitions);


for iT = 1:length(TAGS)
     pos = findcellstr(myXLabels(:,1),TAGS{iT});
     if pos == 0
          xvalues(iT) = iT;
          xlabels{iT} = num2str(iT);
      else
          xvalues(iT)  = cell2mat(myXLabels(pos,3));
          xlabels{iT} = myXLabels(pos,2);
      end
end

% if iscell(xlabels)
%     xlabels = [xlabels{:}];
% end

