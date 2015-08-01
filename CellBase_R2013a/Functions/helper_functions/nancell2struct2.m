function s = nancell2struct2(s)
%NANCELL2STRUCT2   Convert cell of structures with NaNs to struct.
%   S = NANCELL2STRUCT2(C) converts a cell array of structures with the
%   same fields (C) into a struct (S). Fields are filled with NaNs for NaNs
%   in C.
%
%   See also CELL2STRUCT and NANCELL2MAT.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   22-July-2013

%   Edit log: BH 7/22/13

% Convert cell of structures with NaNs to struct
naninx = cellfun(@(s)~isstruct(s),s);   % indices for NaNs
ginx = find(~naninx,1,'first');
fld = fieldnames(s{ginx});   % fieldnames
str = [fld'; repmat({''',NaN,'''},size(fld'))];
str = str(:)';
str = cell2mat(str);
str = ['''' str(1:end-2)];
str = ['struct(' str ');'];
s(naninx) = {eval(str)};   % replace NaNs with empty structures
s = [s{:}];   % convert to struct