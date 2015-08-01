function [D, P, SE] = rocarea(x,y,varargin)
%ROCAREA   Receiver-operator characteristic.
%   D = ROCAREA(X,Y) calculates area under ROC curve (D) for the samples X
%   and Y.
%
%   [D P SE] = ROCAREA(X,Y,'BOOTSTRAP',N) calculates permutation test
%   p-value and bootstrap standard error (N, number of bootstrap samples).
%
%   Optional input parameter-value paits (with default values):
%       'bootstrap', 0 - size of bootstrap sample; 0 corresponds to no
%           bootstrap analysis
%       'transform', 'none' - 'swap': rescales results between 0.5 - 1
%           'scale' - rescales results between -1 to 1
%       'display', false - controls plotting of the ROC curve

%   Edit log: AK 2/2002; AK 4/2005, BH 10/29/13

% Default arguments
prs = inputParser;
addRequired(prs,'x',@isnumeric)   % first sample
addRequired(prs,'y',@isnumeric)   % second sample
addParamValue(prs,'bootstrap',0,@isnumeric)   % size of bootstrap sample
addParamValue(prs,'transform','none',@(s)ismember(s,{'none' 'swap' 'scale'}))   % rescaling
addParamValue(prs,'display',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
parse(prs,x,y,varargin{:})
g = prs.Results;

% Binning
x = x(:);   % convert to column vector
y = y(:);
nbin = ceil(max(length(x)*1.2,length(y)*1.2));   % bin number
mn = min([x; y]);
mx = max([x; y]);
if isequal(mn,mx);   % all numbers in x and y are the same: no discrimination
    D = 0;
    P = 1;
    SE = 0;
    return
end
bin_size = (mx - mn) / nbin;   % bin size
bins = mn-bin_size:bin_size:mx+bin_size;
Lx = length(x);
Ly = length(y);

% AUROC
[D cdfx cdfy] = auc(x,y,Lx,Ly,bins);

% Bootstrap
if g.bootstrap > 0
    z = [x; y];
    Dboot = nan(1,g.bootstrap);
    for k = 1:g.bootstrap
        order = round(rand(1,Lx+Ly)*(Lx+Ly-1))+1;  % resample
        px = z(order(1:Lx));
        py = z(order(Lx+1:Lx+Ly));
        Dboot(k)= auc(px,py,Lx,Ly,bins);   % bootstrap ROC
    end
    
    % p-value
    P = iprctile(Dboot,D);
    if D > mean(Dboot)   % decide which side it should be on
        P = 1 - P;
    end
    
    % SE
    SE = std(Dboot);   % bootstrap standard error
end

% Rescale
switch g.transform
    case 'swap'
        D = abs(D-0.5) + 0.5;   % 'swap'
    case 'scale'
        D = 2 * (D - 0.5);   % 'scale'
end

% Plot
if g.display
    figure
    hold on
    plot(cdfx,cdfy,'b','LineWidth',2);
    plot([0 1],[0 1],'k');
    xlabel('x');
    ylabel('y');
    title(num2str(D));
end

% -------------------------------------------------------------------------
function [D cdfx cdfy] = auc(x,y,Lx,Ly,bins)

% Distributions
px = histc(x,bins);   % distribution of first sample
px = px(1:end-1);
py = histc(y,bins);   % distribution of second sample
py = py(1:end-1);
cdfx = cumsum(px) / Lx;   % CDF of first sample
cdfy = cumsum(py) / Ly;   % CDF of second sample

% AUROC
if isempty(cdfx) || isempty(cdfy)
    D = NaN;
else
    D = trapz(cdfx,cdfy);   % ROC area under the curve
end