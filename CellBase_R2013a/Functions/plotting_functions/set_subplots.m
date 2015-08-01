function H = set_subplots(nrows,ncols,interw,interh,varargin)
%SET_SUBPLOTS   Make subplots with defined distances.
%	H = SET_SUBPLOTS(NROWS,NCOLS,INTERW,INTERH,PARAM,VALUE) creates
%	subplots in an NROWS x NCOLS arrangement with horizontal and vertical
%	distances between the plots given in INTERW and INTERH, respectively.
%	The parameter, value pairs can be used to specify axes properties
%	similar to SUBPLOT.
%
%   See also SUBPLOT.

%   Edit log: BH 5/7/12

% Input arguments
if nargin < 1
    ncols = 4;  % number of columns
end
if nargin < 2
    nrows = 5;  % number of rows
end
if nargin < 3
    interw = 0.02;  % distance between columns
elseif isempty(interw)
    interw = 0.02;
end
if nargin < 4
    interh = 0.02;  % distance between rows
elseif isempty(interh)
    interh = 0.02;
end

% Subplot width and hight
pw = (1 - ((ncols + 1) * interw)) / ncols;  % width of subplots
ph = (1 - ((nrows + 1) * interh)) / nrows;  % hight of subplots

% Create subplots
H = nan(1,nrows*ncols);
try
    for row = 1:nrows
        for col = 1:ncols
            co = (row - 1) * ncols + col;
            left = col * interw + (col - 1) * pw;
            bottom = 1 - (row * interh + row * ph);
            H(co) = subplot(nrows,ncols,co,'Position',[left bottom pw ph],varargin{:});  % make subplot
        end
    end
catch ME    % if positioning fails, make subplots the usual way
    warning(ME.message)
    disp('Calling ''subplot'' without position argument.')
    for row = 1:nrows
        for col = 1:ncols
            co = (row - 1) * ncols + col;
            H(co) = subplot(nrows,ncols,co,varargin{:});
        end
    end
end