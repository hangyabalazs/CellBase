function ind = match_list(TOMATCH, LIST)
%MATCH_LIST   Find exact matches of a string in a cell array of strings.
%   I = MATCH_LIST(PATTERNS,STRINGARRAY) finds exact matches of PATTERNS in
%   STRINGARRAY.
%
%   See also STRMATCH.

%   Edit log: AK 7/1,  BH 1/28/13

% Find matches
NumPatterns = length(TOMATCH);   % number of patterns to match
ind = nan(1,NumPatterns);   % preallocate output
for i = 1:length(TOMATCH)
    ind(i) = strmatch(TOMATCH{i},LIST,'exact');   % find exact matches
end