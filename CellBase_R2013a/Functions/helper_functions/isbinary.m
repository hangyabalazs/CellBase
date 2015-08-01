function TrueFalse = isbinary(x)
%ISBINARY   Determine whether an array is binary.
%   TF = ISBINARY(X) returns 1 if X consists of 0 and 1 values, and 0 
%   otherwise.
%
%   See also DEC2BIN and BIN2DEC.

%   Edit log: BH 1/28/13

% Determine whether the input is binary
TrueFalse =  isempty(setdiff(unique(x),[0 1]));