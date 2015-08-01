function X_xls = formatforExcel(X)
%FORMATFOREXCEL   Helper function for writing numbers to Excel files.
%   X_XLS = FORMATFOREXCEL(X) changes special numbers as Inf, -Inf and NaN
%   to strings. This is useful when numbers are written to Excel to
%   preserve special numbers.
%
%   See also XLSWRITE.

% Convert special numbers to string
if numel(X) == 1   % scalar input
    if ismember(X,[Inf,-Inf]) || isnan(X)      % Inf is converted to 65535 in Excel
        X_xls = {num2str(X)};  % output is converted to cell
    else
        X_xls = X;  % output = input
    end
else   % matrix input
    infinx = cellfun(@(s)~isempty(s)&&isnumeric(s)&&ismember(s(1),[Inf,-Inf]),X(:));   % find Inf and -Inf
    naninx = cellfun(@(s)~isempty(s)&&isnumeric(s)&&isnan(s(1)),X(:));   % find NaNs
    if any(infinx) || any(naninx)
        X_xls = X;
        X_xls(infinx|naninx) = cellfun(@(s){num2str(s)},X(infinx|naninx));
    else
        X_xls = X;  % output = input
    end
end