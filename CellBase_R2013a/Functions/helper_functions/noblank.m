function s1 = noblank(s)
%NOBLANK   Remove blanks.
%   NOBLANK(S) removes trailing blanks from string S.
%
%   NOBLANK(C), when C is a cell array of strings, removes the trailing
%   blanks from each element of C.
%
%   See also PARTITION_TRIALS.

%   Edit log: AK 07/05, BH 6/23/11

% Remove blanks
if isempty(s)
   s1 = s([]);
else
   if isstr(s) 
      s1 = s(~isspace(s));
   elseif  iscellstr(s)
    for iS = 1:length(s)
         j = s{iS};
         s1{iS} = j(~isspace(j));
    end
   else
      warning('MATLAB:deblank:NonStringInput','Input must be a string.')
  end
end