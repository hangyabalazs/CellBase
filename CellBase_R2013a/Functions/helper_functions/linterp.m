function yi = linterp(x,y,xi)
%LINTERP   1D linear interpolation.
%   YI = LINTERP(X,Y,XI) interpolates Y = f(XI), where f is given in (X,Y).
%
%   See also INTERP1Q.

% Interpolation
lx = length(xi);
lex = length(x);
yi = zeros(1,lx);
for k = 1:lx    % linear interpolation
    inx1 = find(x<=xi(k),1,'last');
    inx1(inx1==lex) = inx1(inx1==lex) - 1;
    inx2 = inx1 + 1;
    yi(k) = y(inx1) + (y(inx2) - y(inx1)) * (xi(k) - x(inx1)) / (x(inx2) - x(inx1));
end