function varargout = errorshade(X,Y,SE,varargin)
%ERRORSHADE   Shaded error bar plot.
%   H = ERRORSHADE(X,Y,SD,LINECOLOR) plots Y against X in LineColor. A
%   shaded error bar is added between Y-SD and Y+SD. Output returns the
%   handles for the line and the shaded error bar. If SD is a 2-by-N
%   matrix, the two rows are used as upper and lower bound of the shade.
%
%   H = ERRORSHADE(X,Y,SD,'LINECOLOR',LCLR,'LINEWIDTH',LW,'SHADECOLOR',SCLR)
%   - optional arguments LineWidth and ShadeColor default to 2 and grey.
%   LineColor defaults to blue - you can use shortend styles as well ('r:').
%
%   See also PLOT_TIMECOURSE.

%   Edit log: AK 6/2004, BH 6/27/11, BH 1/31/12, BH 3/11/12

% Assign defaults if no arguments passed
default_args = { ...
        'LineWidth'         2; ...
        'LineColor'         'b'; ...
        'LineStyle'         '-'; ...
        'ShadeColor'        'gray'; ...
        'FaceAlpha'         0.3;  ...
        'PlotDashedTime'    NaN;  ...
        'DashedLineStyle'   '--'; ...
    };
if isequal(length(varargin),1)
    varargin{2} = varargin{1};
    varargin{1} = 'LineColor';
end
[g, error] = parse_args(default_args,varargin{:});
if ~isnumeric(g.ShadeColor)
    switch lower(g.ShadeColor)
        case 'same'
            g.ShadeColor = min(g.LineColor(:)' + [0.4 0.4 0.4],1);
        case {'gray','grey'}
            g.ShadeColor = [0.6 0.6 0.6];
    end
end

% Input argument check
if (min(size(SE)) > 1) 
   E = SE;    
else
   E(1,:) = Y - SE;
   E(2,:) = Y + SE;
end 
posNoNaN = find(~isnan(Y));
X = X(posNoNaN);
Y = Y(posNoNaN);
E = E(:,posNoNaN);

% Plot
hold on;
p1 = patch([X X(end:-1:1)],[E(1,:) E(2,end:-1:1)],g.ShadeColor);
hAnnotation = get(p1,'Annotation');
hLegendEntry = get(hAnnotation','LegendInformation');
set(hLegendEntry,'IconDisplayStyle','off')
set(p1,'LineStyle','none','EdgeColor',g.ShadeColor,'FaceColor',g.ShadeColor,...
    'FaceAlpha',g.FaceAlpha);
if length(g.LineColor) > 1 && ~isnumeric(g.LineColor)
    p2 = plot(X,Y,'Color',g.LineColor(1),'LineWidth',g.LineWidth);
    set(p2,'linestyle',g.LineColor(2:end))
else
    if ~isnan(g.PlotDashedTime)
        posDASHED = nearest(X,g.PlotDashedTime);
        p2(1) = plot(X(1:posDASHED),Y(1:posDASHED),'Color',g.LineColor,'LineWidth',g.LineWidth,'LineStyle',g.LineStyle);
        p2(2) = plot(X(posDASHED:end),Y(posDASHED:end),'Color',g.LineColor,'LineWidth',g.LineWidth,'LineStyle',g.DashedLineStyle);
    else        
        p2 = plot(X,Y,'Color',g.LineColor,'LineWidth',g.LineWidth,'LineStyle',g.LineStyle);
    end
end

% Output
if nargout == 1 
    varargout{1} = [p1 p2];
end