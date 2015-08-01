function s = nancell2mat(s)
%NANCELL2MAT   Convert cell of equal size matrices with NaNs to double
%   S = NANCELL2MAT(C) converts a cell array of matrices with equal size 
%   (C) into a double (S). NaNs in C are padded to the proper size.
%
%   See also CELL2MAT and NANCELL2STRUCT.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   2-Oct-2013

%   Edit log: BH 8/2/13

% Convert cell of structures with NaNs to struct
naninx = cellfun(@(s)all(isnan(s))&isequal(numel(s),1),s);   % indices for NaNs
ginx = find(~naninx,1,'first');
sz = size(s{ginx});   % size of the matrices
s(naninx) = {nan(sz)};   % make NaN matrices of the same size
s = cell2mat(s);   % convert to double