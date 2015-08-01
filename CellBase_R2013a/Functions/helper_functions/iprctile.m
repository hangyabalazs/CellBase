function p = iprctile(x,value)
%IPRCTILE   Percentile.
%   P = IPRCTILE(X,PCT) returns the PCT (between 0 and 1) percentile of X.
%
%   See also PRCTILE.

%   Edit log: AK 2/2005, BH 1/28/13

% Convert to column vector
x = x(:)';

% If all values are the same, return that value
if length(unique(x)) == 1
    p = (unique(x) == value);
    return
end

% Get percentile by intepolation
xx = sort(x);   % sorted input
Lx = length(x);   % length of input
ind = find(diff(xx)~=0);  % indices where the input changes
ind = [ind+1 ind(end)+2];
xx = [-Inf xx Inf];
q = [0 ((0.5:Lx-0.5)/Lx) 1];  % predefined percentiles for interpolation
try
    p = interp1(xx(ind),q(ind),value,'linear','extrap');   % interpolate percentile
    p = max(0,min(1,p));
catch ME
    disp(ME.message)   % interpolation failed
    warning('ipercentile:interpolationProblem','Error in IPERCENTILE.')
    p = NaN;
end