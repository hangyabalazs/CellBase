function sX = standardize(X)
%STANDARDIZE   Standardization.
%   SX = STANDARDIZE(X) standardizes X to zero mean and unitary standard
%   deviation.
%
%   See also ZSCORE.

sX = (X - nanmean(X)) / nanstd(X);