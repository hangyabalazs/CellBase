function hout = fstamp(str,varargin)
%FSTAMP   Put a title below all subplots.
%   FSTAMP('TEXT',FONTSIZE) adds text to the bottom of the figure below all
%   subplots. Use it for stamping filenames. 
%
%   It uses an optional FONTSIZE argument. If you pass fontsize multiplied
%   by 10, you get bold fonts.

%   Edit log: Drea Thomas 6/15/95 drea@mathworks.com, AK 3/22/01, BH
%   6/27/11

% Amount of the figure window devoted to subplots
plotregion = 0.91;  
str = strrep(str,'_','-');

% X & Y position of title in normalized coordinates
titlexpos  = .92;
titleypos  = .03;

% Fontsize for supertitle
ITALIC = 0;
if (nargin > 1)
   fs = varargin{1};
   if (fs>40)
       ITALIC = 1;
       fs = fs / 10;
   end
else
   fs = get(gcf,'defaultaxesfontsize') + 2;
end

if nargin > 2
    titlepos = varargin{2};
else
    titlepos = 'bottom_right';
end

% Fudge factor to adjust y spacing between subplots
fudge = 0.7;
haold = gca;
figunits = get(gcf,'units');

% Get the (approximate) difference between full height (plot + title
% + xlabel) and bounding rectangle.
if ~strcmp(figunits,'pixels')
    set(gcf,'units','pixels');
    pos = get(gcf,'position');
    set(gcf,'units',figunits);
else
    pos = get(gcf,'position');
end
ff = (fs-4) * 1.27 * 5 / pos(4) * fudge;
% The 5 here reflects about 3 characters of height below an axis and 2
% above. 1.27 is pixels per point.

% Determine the bounding rectange for all the plots
h = findobj(gcf,'Type','axes');  % Change suggested by Stacy J. Hills
max_y = 0;
min_y = 1;

oldtitle = 0;
for i = 1:length(h)
	if ~strcmp(get(h(i),'Tag'),'fstamp')
		pos = get(h(i),'pos');
    else
		oldtitle = h(i);
	end
end

if max_y > plotregion
	scale = (plotregion-min_y) / (max_y-min_y);
	for i = 1:length(h)
		pos = get(h(i),'position');
		pos(2) = (pos(2) - min_y) * scale + min_y;
		pos(4) = pos(4) * scale - (1 - scale) * ff / 5 * 3;
		set(h(i),'position',pos);
	end
end

np = get(gcf,'nextplot');
set(gcf,'nextplot','add');
if (oldtitle),
	delete(oldtitle);
end
ha = axes('pos',[0 1 1 1],'visible','off','Tag','fstamp');
ht = text(titlexpos,titleypos-1,str);

% Positioning
switch titlepos
    case 'bottom-right'
        set(ht,'horizontalalignment','right','fontsize',fs);
    case 'top-center'
        titlexpos = 0.48;
        titleypos = -0.02;
        set(ht,'Position',[titlexpos titleypos],'VerticalAlignment','top','HorizontalAlignment','center')
end
if (ITALIC)
  set(ht,'FontAngle','italic');
end
set(gcf,'nextplot',np);
axes(haold);

% Output
if nargout
	hout = ht;
end