function D = KLdist(P,Q)
%KLDIST   Kullbach-Leibler distance.
%   D = KLDIST(P,Q) calculates the Kullbach-Leibler distance (information
%   divergence) of the two input distributions.
%
%   See also JSDIV, FDIV, HDISC, BDIST, CHISQUAREDIV, VARDIST and 
%   HARMONICMEAN.

% Input argument check
error(nargchk(2,2,nargin))
if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
    error('Input arguments must be probability distributions.')
end
if ~isequal(size(P),size(Q))
    error('Input distributions must be of the same size.')
end

% KL-distance
P2 = P(P.*Q>0);     % restrict to the common support
Q2 = Q(P.*Q>0);
P2 = P2 / sum(P2);  % renormalize
Q2 = Q2 / sum(Q2);

D = sum(P2.*log(P2./Q2));

% Alternative way of computation:
% HPQ = -sum(P2.*log(Q2));      % cross-entropy
% HP = -sum(P2.*log(P2));       % entropy
% D = HPQ - HP;