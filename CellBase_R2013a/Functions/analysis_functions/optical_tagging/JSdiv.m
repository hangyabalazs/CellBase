function D = JSdiv(P,Q)
%JSDIV   Jensen-Shannon divergence.
%   D = JSDIV(P,Q) calculates the Jensen-Shannon divergence of the two 
%   input distributions.
%
%   See also KLDIST, FDIV, HDISC, BDIST, CHISQUAREDIV, VARDIST and 
%   HARMONICMEAN.

% Input argument check
error(nargchk(2,2,nargin))
if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
    error('Input arguments must be probability distributions.')
end
if ~isequal(size(P),size(Q))
    error('Input distributions must be of the same size.')
end

% JS-divergence
M = (P + Q) / 2;
D1 = KLdist(P,M);
D2 = KLdist(Q,M);
D = (D1 + D2) / 2;