function [p, h, stats] = b_signrank2(x,y,varargin)
% one-sided

% Check most of the inputs now
alpha = 0.05;
if nargin>2 && isnumeric(varargin{1})
   % Grandfathered syntax:  signrank(x,y,alpha)
   alpha = varargin{1};
   varargin(1) = [];
end
oknames = {'alpha' 'method'};
dflts   = {alpha   ''};
[eid,emsg,alpha,method] = b_statgetargs(oknames,dflts,varargin{:});
if ~isempty(eid)
   error(sprintf('stats:signrank:%s',eid),emsg);
end

if ~isscalar(alpha)
   error('stats:signrank:BadAlpha','SIGNRANK requires a scalar ALPHA value.');
end
if ~isnumeric(alpha) || isnan(alpha) || (alpha <= 0) || (alpha >= 1)
   error('stats:signrank:BadAlpha','SIGNRANK requires 0 < ALPHA < 1.');
end

if nargin < 2 || isempty(y)
    y = zeros(size(x));
elseif isscalar(y)
    y = repmat(y, size(x));
end

if ~isvector(x) || ~isvector(y)
    error('stats:signrank:InvalidData',...
          'SIGNRANK requires vector rather than matrix data.');
elseif numel(x) ~= numel(y)
    error('stats:signrank:InputSizeMismatch',...
    'SIGNRANK requires the data vectors to have the same number of elements.');
end

diffxy = x(:) - y(:);

% Remove missing data
diffxy(isnan(diffxy)) = [];
if (length(diffxy)==0)
   error('stats:signrank:NotEnoughData','No data remaining after removal of NaNs.');
end

nodiff = find(diffxy == 0);
diffxy(nodiff) = [];
n = length(diffxy);

if (n == 0)         % degenerate case, all ties
    p = 1;
    if (nargout > 1)
        h = 0;
        if (nargout > 2)
            stats.signedrank = 0;
        end
    end
    return
end

% Now deal with the method argument
if isempty(method)
   if n<=15
      method = 'exact';
   else
      method = 'approximate';
   end
elseif ischar(method)
   okmethods = {'exact' 'approximate' 'oldexact'};
   j = strmatch(lower(method),okmethods);
   if isempty(j)
      error('stats:signrank:BadMethod',...
            'METHOD must be ''exact'' or ''approximate''.');
   end
   method = okmethods{j};
else
   error('stats:signrank:BadMethod',...
         'METHOD must be ''exact'' or ''approximate''.');
end

% Find negative differences and ranks of absolute differences
neg = find(diffxy<0);
[tierank, tieadj] = tiedrank(abs(diffxy));

% Compute signed rank statistic (most extreme version)
w = sum(tierank(neg));
w = min(w, n*(n+1)/2-w);

if isequal(method,'approximate')
    z = (w-n*(n+1)/4) / sqrt((n*(n+1)*(2*n+1) - tieadj)/24);
    p = normcdf(z,0,1);
    if (nargout > 2)
        stats.zval = z;
    end
elseif isequal(method,'oldexact')
    % Enumerates all possibilities and does not adjust for ties
    allposs = (ff2n(n))';
    idx = (1:n)';
    idx = idx(:,ones(2.^n,1));
    pranks = sum(allposs.*idx,1);
    tail = length(find(pranks <= w)); % one side.

    % Avoid p>1 if w is in the middle and is double-counted
    p = min(1, tail./(2.^n));
else % isequal(method,'exact')
    p = b_statsrexact(tierank,w);
    p = min(1, p);   % one-sided, don't double-count the middle value
end

if nargout > 1
    h = (p<=alpha);
    if (nargout > 2)
        stats.signedrank = w;
    end
end
