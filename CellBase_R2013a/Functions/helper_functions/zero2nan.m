function Y = zero2nan(X)
%ZERO2NAN   Replace zeros with NaNs.
%   Y = ZERO2NAN(X) replaces zeros in X with NaNs.
%
%   See also NAN2ZERO.

Y = X;
Y(X==0) = NaN;