function  yn = cellcmp(cell1, cell2)
%CELLCMP   Compare two cell arrays.
%   YN = CELLCMP(CELL1,CELL2) compares two cell arrays and returns 1 if the
%   cell arrays are the same and 0 otherwise. It works only on cell arrays
%   of matrices.
%
%   See also ISEQUAL.

%   Edit log: BH 1/28/13

% Compare two cell arrays
if iscell(cell1) && iscell(cell2)  % only if inputs are cells
    L1 = length(cell1);
    L2 = length(cell2);
    if L1 == L2   % compare length
        yn = 1;
        for i = 1:L1   % compare each cell one-by-one
            c1 = cell2mat(cell1(i));   % convert to double
            c2 = cell2mat(cell2(i));
            if length(c1) ~= length(c2)   % compare length
                yn = 0;
            else
                yn = yn * prod(double(c1==c2));   % compare content
            end
        end
        return
    else
        yn = 0;
        return
    end
else
    yn = 0;
end