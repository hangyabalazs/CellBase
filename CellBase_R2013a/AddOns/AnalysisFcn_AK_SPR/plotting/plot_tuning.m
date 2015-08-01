function vargout = plot_tuning(X,Y,params,varargin)
% 

% define_fun = @defineLabelsColors;

default_args = params;
default_args.TitleStr = '';
default_args.XLabel = 'X';
default_args.YLabel = 'Y';
default_args.Legend = '';
default_args.Color  = 'b';
default_args.MarkerSize = 8;
default_args.LineWidth = 3;
default_args.LineColor = 'k';
default_args.Line = 'on';
default_args.XTickLabels = '';
default_args.Limits = [];
[par, error] = parse_args(default_args,varargin{:});

%-------------------------------------------------------------
NumPlots = [];

CROSSBAR = 0;
hebar=1; heline=2;
XAUTO = 0; 
if isnan(X) | isempty(X)
    X = [1:length(Y)]';
    XAUTO = 1;
%elseif size(X,2) == 2
elseif size(X,1) == 2
    CROSSBAR=1;
    hebar=[1 2]; heline=3;
end
if min(size(X))==1
    X=X(:);  
end

if sum(isnan(Y(:,1))) > 0
    %there are NaNs
    pos = find(~isnan(Y(:,1)));
    Y = Y(pos,:);
    if size(X,2) == 2
        X = X(pos,:);
    else
        X = X(pos);
    end
    if XAUTO
        X = [1:length(Y)]';
    end
    if ~isempty(par.XTickLabels)
        par.XTickLabels = par.XTickLabels(pos);
    end
end

if length(par.LineWidth)==1
    par.LineWidth(2) = par.LineWidth(1);
end

f0=focusfigure(par.FigureNum);
hold on;
if iscell(par.Color)  
    for iX=1:length(X(:,1))
        if CROSSBAR
            he(:,iX) = errorbar2(X(iX,1),X(iX,2),Y(iX,1),Y(iX,2),'t0');            
        elseif size(Y,2) == 2
            he(:,iX) = errorbar(X(iX,1),Y(iX,1),Y(iX,2),'t0');            
        else
            he(:,iX) = plot(X(iX,1),Y(iX,1),'-');      
        end
        
        if par.MarkerSize>0
            hm(iX) = plot(X(iX,1),Y(iX,1),'o');
            set(hm(iX),'MarkerFaceColor',par.Color{iX},'MarkerEdgeColor','none','MarkerSize',par.MarkerSize);
        end
        set(he(hebar,iX),'Color',par.Color{iX});
    end
    set(he(hebar),'LineWidth',par.LineWidth(2));
    delete(he(heline));
else
    if CROSSBAR
        he = errorbar2(X(:,1),X(:,2),Y(:,1),Y(:,2),'t0');      
    else
        he = errorbar(X(:,1),Y(:,1),Y(:,2),'t0');
    end
    if par.MarkerSize>0
        hm = plot(X(:,1),Y(:,1),'o');
        set(hm,'MarkerFaceColor',par.Color,'MarkerEdgeColor','none','MarkerSize',par.MarkerSize);
    end     
    set(he,'Color',par.Color);
    set(he(hebar),'LineWidth',par.LineWidth(2));
    delete(he(heline));
end

if strcmpi(par.Line,'on')
    pos=~isnan(Y(:,1));
    hl = plot(X(pos,1),Y(pos,1),'Color',par.LineColor);
    set(hl,'LineWidth',par.LineWidth(1));
    if ~iscell(par.Color) 
        set(hl,'Color',par.Color);
    end
end


min_x = min(X(:,1)); max_x = max(X(:,1));
min_y = min(Y(:,1)); max_y = max(Y(:,1));
margin_x = abs(diff(X([1 end],1)))/5;
margin_y = (max_y-min_y)/5;
axis([min_x-margin_x max_x+margin_x+2*eps min_y-margin_y max_y+margin_y+2*eps]);

if (length(unique(diff(X(:,1)))) == 1) & ~isempty(find(X))
    if diff(X(:,1)) < 0
        'oops - plot_tuning line 108'
    else
        set(gca,'XTick',X(:,1));
    end
end

if ~isempty(par.XTickLabels)
    set(gca,'XTickLabel',par.XTickLabels);
end
if ~isempty(par.XTickLabels)
    set(gca,'XTickLabel',par.XTickLabels);
end

%-------------------------------------------------------
l(1)=ylabel(par.YLabel);
l(2)=xlabel(par.XLabel);
li=2;
if ~isempty(par.TitleStr)
    li=li+1;
    l(li)=title(par.TitleStr);
end

if ~isempty(par.Legend)
    l_leg=legend(par.Legend(NumPlots),1);
    set(l_leg,'box','off','FontSize',8,'Color','none','XColor','w','YColor','w');
end
%set(gca,'XTickLabelMode','auto','XColor','k','YTickLabelMode','auto','YColor','k','YTickMode','auto');
%
alim = axis;
if ~isempty(par.Limits)
    pos=~isnan(par.Limits);
    alim(pos) = par.Limits(pos);   
end

if  size(X,2) == 2
    xmin=min([X(:,1)-X(:,2); X(:,1)+X(:,2)]);
    xmax=max([X(:,1)-X(:,2); X(:,1)+X(:,2)]);
else
    xmin=min(X(:));
    xmax=max(X(:));
end
if  size(Y,2) == 2
    ymin=min([Y(:,1)-Y(:,2); Y(:,1)+Y(:,2)]);
    ymax=max([Y(:,1)-Y(:,2); Y(:,1)+Y(:,2)]);
else
    ymin=min(Y(:));
    ymax=max(Y(:));
end
dx=(xmax-xmin)/8;
xmin=xmin-dx; xmax=xmax+dx;
%

dy=(ymax-ymin)/8;
if dy==0
    dy=2*eps;
end
ymin=ymin-dy; ymax=ymax+dy;
alim(1) = min(alim(1),xmin);
alim(2) = max(alim(2),xmax);
alim(3) = min(alim(3),ymin);
alim(4) = max(alim(4),ymax);

axis(alim);

setmyplot(gca,[l]);
if nargout > 0
    vargout = f0;
end

