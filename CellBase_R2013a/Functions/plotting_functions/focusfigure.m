function fhandle = focusfigure(H,varargin)
%FOCUSFIGURE   Bring figure and subplot into focus.
%   FHANDLE = FOCUSFIGURE(H) sets current figure or subplot to H. 
%
%   FOCUSFIGURE(H,'CREATE') creates new figure.
%
%   See also PLOT_TIMECOURSE.

%   Edit log: AK 2/07, BH 6/27/11

% Input argument check
CREATE = 0;

if (min(ishandle(H)) < 1) % not handle    
    if isnumeric(H) && (nargin > 1) && strcmpi(varargin{1},'create')
        CREATE = 1;
    else
        error('FOCUSFIGURE: Invalid handle passed.')
    end
end

% Set current figure or subplot
if CREATE   % user asked to create figure
    fhandle = figure(H);
    clf;
else
    lenhan = length(H);
    fhandle = zeros(1,lenhan);
    for iH = 1:lenhan
        if strcmp(get(H(iH),'Type'),'figure')
            figure(H(iH));
            fhandle(iH) = H(iH);
        elseif strcmp(get(H(iH),'Type'),'axes')
            subplot(H(iH));
            fhandle(iH) = H(iH);
        end
    end    % iH
end