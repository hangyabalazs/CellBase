function [X2 S] = smooth(X,str,wn)
%SMOOTH   Smooth with a moving average.
%   [O S] = SMOOTH(X,M,WN) performs smoothing on X by averaging in a
%   sliding window of WN size. M describes the data type, i.e. 'linear' or 
%   'circular'. Output is returned in O with standard error of the
%   smoothing in S.
%
%   See also CONV.

n = wn;
nn = (n - 1) / 2;
m = length(X);
X2 = zeros(size(X));
S = zeros(size(X));
switch str
    case 'linear'
        for k = 1:m
            ind1 = max(1,k-nn);
            ind2 = min(k+nn,m);
            X2(k) = mean(X(ind1:ind2));
            S(k) = std(X(ind1:ind2));
        end
    case 'circular'
        for k = 1:m
            if k - nn < 1
                X2(k) = mean([X(mod2(k-nn,m):m); X(1:k+nn)]);
                S(k) = std([X(mod2(k-nn,m):m); X(1:k+nn)]) / sqrt(n);
            elseif k + nn > m
                X2(k) = mean([X(k-nn:m); X(1:mod2(k+nn,m))]);
                S(k) = std([X(k-nn:m); X(1:mod2(k+nn,m))]) / sqrt(n);
            else
                X2(k) = mean(X(k-nn:k+nn));
                S(k) = std(X(k-nn:k+nn)) / sqrt(n);
            end
        end
end