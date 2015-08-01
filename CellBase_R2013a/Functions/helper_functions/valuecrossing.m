function t = valuecrossing(x,y,v,opt)
%VALUECROSSING   Interpolate value-crossings of a function.
%   T = VALUECROSSING(X,Y,V) returns the V-crossings of the function (X,Y).
%   VALUECROSSING uses linear interpolation.
%   T = VALUECROSSING(X,Y,V,OPT) determines upwards, downwards or all
%   value-crossings using 'up', 'down' or 'both' options, respectively.

%   Balazs Hangya
%   Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   New York 11724, USA
%   balazs.cshl@gmail.com
%
%   Institute of Experimental Medicine
%   Szigony street 43, Budapest 1083, Hungary
%   hangyab@koki.hu

% Input argument check
error(nargchk(3,4,nargin))
if nargin == 3
    opt = 'both';
end

% Interpolate crossings
lup = y <= v & [y(2:end) v] > v;
ldown = y >= v & [y(2:end) v] < v;
switch opt
    case 'up'
        fnd = find(lup);
    case 'down'
        fnd = find(ldown);
    case 'both'
        fnd = find(lup|ldown);
end
xf = x(fnd);
xff = x(fnd+1);
yf = y(fnd);
yff = y(fnd+1);
mpl = (v - yf) ./ (yff - yf);
t = xf + (xff - xf) .* mpl;