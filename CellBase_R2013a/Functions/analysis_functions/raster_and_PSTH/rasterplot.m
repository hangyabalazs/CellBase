function varargout = rasterplot(varargin)
%RASTERPLOT   Raster plot.
%   H = RASTERPLOT(BINRASTER) plots spike raster from BINRASTER. It shows a
%   'histogram' view if spikes are overlapping on the plot (more spikes
%   appear in darker colors). Zoom is implemented as buttondown function
%   (usage is the same as the usual zoom).
%
%   H = RASTERPLOT(BINRASTER,TIME) uses the TIME array for specifying
%   x-coordinates.
%
%   RASTERPLOT(BINRASTER,TIME,H) and RASTERPLOT(BINRASTER,[],H) plots in
%   the figure with handle H.
%
%   See also STIMES2BINRASTER, VIEWCELL2B and PLOT_RASTER2A.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   09-Oct-2011

%   Edit log: BH, 9/10/2011

% Launch
if ~ischar(varargin{1})
    
    binraster = varargin{1};
    
    if length(varargin) < 3
        H = figure;
    else
        H = varargin{3};
        if ~ishandle(H) || (~isequal(get(H,'Type'),'axes') && ~isequal(get(H,'Type'),'figure'))
            error('rasterplot:InvHandle','Invalid handle. Third input argument should be a valid figure or axes handle.')
        end
        if isequal(get(H,'Type'),'axes')
            H = ancestor(H,'figure');
        end
    end
    
    [varargout{1:nargout}] = openraster(H,binraster,varargin{2:end});
else    % invoke named subfunction or callback
	[varargout{1:nargout}] = feval(varargin{:});
end

% -------------------------------------------------------------------------
function H = openraster(H,binraster,varargin)

% Prepare figure
ex = false;
if length(varargin) > 1
    ex = true;   % if the figure existed, ResizeFcn is not called and 'main' should be called directly
    if isequal(get(varargin{2},'Type'),'axes')
        A = varargin{2};
        axes(A);   %#ok<MAXES>
    else
        figure(H);
        A = axes;
    end
else
    A = axes;
end
cmp = colormap('bone');
colormap(flipud(cmp))
set(A,'Units','pixels')
set(H,'Units','pixels')
pos = get(A,'Position');

% Position
lft = pos(1);
btm = pos(2);
wdth = floor(pos(3));
hght = floor(pos(4));
setappdata(H,'left',lft);
setappdata(H,'bottom',btm);
setappdata(H,'width',wdth)
setappdata(H,'hight',hght)

% Data
szb = size(binraster);
time = 1:szb(2);
if length(varargin) >= 1 && ~isempty(varargin{1})   % 'time' input argument
    rtime = varargin{1};
    if ~isequal(length(rtime),length(time))
        error('rasterplot:inputargDimMismatch',...
            'Length of ''time'' input argument should match binraster dimensions.')
    end
else
    rtime = time;
end
tno = 1:szb(1);
setappdata(H,'binraster',binraster)
setappdata(H,'time',time)
setappdata(H,'time0',time)
setappdata(H,'realtime',rtime)
setappdata(H,'realtime0',rtime)
setappdata(H,'tno',tno)
setappdata(H,'tno0',tno)

% Callbacks
rsf = 'rasterplot(''figure_ResizeFcn'',gcbo)';
set(H,'ResizeFcn',rsf)
bdf = 'rasterplot(''figure_ButtonDownFcn'',gcf)';
setappdata(H,'ButtonDownFcn',bdf);

% If the figure existed, ResizeFcn is not called and 'main' should be called directly
if ex
    main(H,binraster);
end

% -------------------------------------------------------------------------
function H = main(H,binraster)

% Position
wdth = getappdata(H,'width');
hght = getappdata(H,'hight');

% Raster plot
szb = size(binraster);
[x y] = find(binraster);
if isempty(x)
    rstr2 = zeros(hght,wdth);
else
    y = ceil(y/szb(2)*wdth);
    if hght > szb(1)    % only compress in y-direction if necessary ('line-view')
        rstr2 = accumarray([x,y],1,[szb(1) wdth]);
    else
        x = ceil(x/szb(1)*hght);
        rstr2 = accumarray([x,y],1,[hght wdth]);
    end
end

% Draw
tno = getappdata(H,'tno');
rtime = getappdata(H,'realtime');
imagesc(rtime,tno,rstr2)

% Set ButtonDownFcn
bdf = getappdata(H,'ButtonDownFcn');
set(gca,'ButtonDownFcn',bdf)
set(allchild(gca),'ButtonDownFcn',bdf)

% Code for 2-point line view
% if 1.5 * szb(1) < hght
%     x = [x; x+1];
%     y = [y; y];
% end
% inx = x > hght;
% x(inx) = [];
% y(inx) = [];

% Code for exact binning, slower
% pxls = combvec(1:hght,1:wdth);
% numcmb = size(pxls,2);
% rstr = nan(hght,wdth);
% for k = 1:numcmb
%     sqr = binraster(floor((pxls(1,k)-1)*szb(1)/hght+1):ceil(pxls(1,k)*szb(1)/hght),...
%         round((pxls(2,k)-1)*szb(2)/wdth+1):round(pxls(2,k)*szb(2)/wdth));
%     rstr(pxls(1,k),pxls(2,k)) = 1 - prod(1-sqr(:));
% end

% Compact code for the previous solution, which is actually even slower
% sqr2 = arrayfun(@(k)binraster(floor((pxls(1,k)-1)*szb(1)/hght+1):ceil(pxls(1,k)*szb(1)/hght),...
%     round((pxls(2,k)-1)*szb(2)/wdth+1):round(pxls(2,k)*szb(2)/wdth)),1:numcmb,'UniformOutput',false);
% prstr = 1 - cellfun(@(cll)prod(1-cll(:)),sqr2);
% rstr2 = reshape(prstr,hght,wdth);

% If we want to start trials from the bottom (zoom should be modified accordingly!)
% imagesc(time,tno,flipud(rstr2))

% -------------------------------------------------------------------------
function figure_ResizeFcn(hObj)     %#ok<DEFNU>

% Set new width and hight, except first open
A = findobj(allchild(hObj),'Type','axes');
if ~isempty(A)
    
    % Get left and bottom position
    lft = getappdata(hObj,'left');
    btm = getappdata(hObj,'bottom');
    
    % Parent can be figure or uipanel
    pta = get(A,'Parent');
    old_units = get(pta,'Units');
    set(pta,'Units','pixels')
    
    % New position
    figpos = get(pta,'Position');
    wdth = floor(figpos(3)-2*lft+25);
    hght = floor(figpos(4)-2*btm+20);
    lbwh = [lft btm wdth hght];
    set(A,'Units','pixels')
    set(A,'Position',lbwh)
    
    % Store new width and hight
    setappdata(hObj,'width',wdth)
    setappdata(hObj,'hight',hght)
    set(pta,'Units',old_units)
end

% Redraw
binraster = getappdata(hObj,'binraster');
main(hObj,binraster);

% -------------------------------------------------------------------------
function figure_ButtonDownFcn(hObj)     %#ok<DEFNU>

% Get variables
A = gca;
binraster = getappdata(hObj,'binraster');
time = getappdata(hObj,'time');
tno = getappdata(hObj,'tno');
tno0 = getappdata(hObj,'tno0');
time0 = getappdata(hObj,'time0');
rtime0 = getappdata(hObj,'realtime0');

% Set axis
seltyp = get(hObj,'SelectionType');
switch seltyp
case 'normal'   % zoom in
    point1 = get(A,'CurrentPoint'); % button down detected
    rbbox;
    point2 = get(A,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    point1(1) = linterp(rtime0,time0,point1(1));
    point2(1) = linterp(rtime0,time0,point2(1));
    if isequal(point1,point2)
        xx = [time(1) time(end)];
        yy = [tno(1) tno(end)];
        xx2 = (abs(xx(2) - xx(1))) / 4;
        yy2 = (abs(yy(2) - yy(1))) / 4;
        xx3(1) = point1(1) - xx2;
        xx3(2) = point1(1) + xx2;
        yy3(1) = point1(2) - yy2;
        yy3(2) = point1(2) + yy2;
        if time0(1) < xx3(1) && time0(end) > xx3(2)
            xlim_data = xx3;
        elseif time0(1) > xx3(1)
            xx_new(1) = time(1);
            xx_new(2) = time(1) + (2 * xx2);
            xlim_data = xx_new;
        elseif time0(end) < xx3(2)
            xx_new(1) = time0(end) - (2 * xx2);
            xx_new(2) = time0(end);
            xlim_data = xx_new;
        end
        if tno0(1) < yy3(1) && tno0(end) > yy3(2)
            ylim_data = yy3;
        elseif tno0(1) > yy3(1)
            yy_new(1) = tno0(1);
            yy_new(2) = tno0(1) + (2 * yy2);
            ylim_data = yy_new;
        elseif tno0(end) < yy3(2)
            yy_new(1) = tno0(end) - (2 * yy2);
            yy_new(2) = tno0(end);
            ylim_data = yy_new;
        end
    else
        xlim_data = [min([point1(1) point2(1)]) max([point1(1) point2(1)])];
        xlim_data(1) = max(xlim_data(1),time0(1));
        xlim_data(2) = min(xlim_data(2),time0(end));
        ylim_data = [min([point1(2) point2(2)]) max([point1(2) point2(2)])];
        ylim_data(1) = max(ylim_data(1),tno0(1));
        ylim_data(2) = min(ylim_data(2),tno0(end));
    end

case 'open'   % set default
    xlim_data = [time0(1) time0(end)];
    ylim_data = [tno0(1) tno0(end)];
    
case 'extend'   % zoom out
    point = get(A,'CurrentPoint'); % button down detected
    point = point(1,1:2);
    point(1) = linterp(rtime0,time0,point(1));
    xx = [time(1) time(end)];
    yy = [tno(1) tno(end)];
    xx2 = abs(xx(2) - xx(1));
    yy2 = abs(yy(2) - yy(1));
    if xx2 > (time0(end) - time0(1)) / 2,
        xx2 = (time0(end) - time0(1)) / 2;
    end
    if yy2 > (abs(tno0(end) - tno0(1))) / 2,
        yy2 = (abs(tno0(end) - tno0(1))) / 2;
    end
    xx3(1) = point(1) - xx2;
    xx3(2) = point(1) + xx2;
    yy3(1) = point(2) - yy2;
    yy3(2) = point(2) + yy2;
    if time0(1) < xx3(1) && time0(end) > xx3(2)
        xlim_data = xx3;
    elseif time0(1) > xx3(1)
        xx_new(1) = time0(1);
        xx_new(2) = time0(1) + (2 * xx2);
        xlim_data = xx_new;
    elseif time0(end) < xx3(2)
        xx_new(1) = time0(end) - (2 * xx2);
        xx_new(2) = time0(end);
        xlim_data = xx_new;
    end
    if tno0(1) < yy3(1) && tno0(end) > yy3(2)
        ylim_data = yy3;
    elseif tno0(1) > yy3(1)
        yy_new(1) = tno0(1);
        yy_new(2) = tno0(1) + (2 * yy2);
        ylim_data = yy_new;
    elseif tno0(end) < yy3(2)
        yy_new(1) = tno0(end) - (2 * yy2);
        yy_new(2) = tno0(end);
        ylim_data = yy_new;
    end
otherwise
    return
end

% Redraw
time2 = round(xlim_data(1)):round(xlim_data(2));
tno2 = round(ylim_data(1)):round(ylim_data(2));
rtime2 = rtime0(time2(1):time2(end));
binraster2 = binraster(tno2,time2);
setappdata(hObj,'realtime',rtime2);
setappdata(hObj,'time',time2);
setappdata(hObj,'tno',tno2);
main(hObj,binraster2);

% -------------------------------------------------------------------------
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