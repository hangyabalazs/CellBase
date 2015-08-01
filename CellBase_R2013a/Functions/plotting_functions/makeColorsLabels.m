function [mylabels, mycolors, mycolors2, mylinestyle] = makeColorsLabels(fun_definitions,TAGS)
%MAKECOLORSLABELS   Colors, labels and line style.
%   [MYLABELS, MYCOLORS, MYCOLORS2, MYLINESTYLE] = MAKECOLORSLABELS(FUN_DEFINITIONS,TAGS)
%   defines colors, labels and line style for plotting in CellBase using
%   the colors defined by FUN_DEFINITIONS.
%
%   See also DEFINELABELSCOLORS_DEFAULT and VIEWCELL2B.

%   Edit log: AK 7/1, BH 6/24/11

% Get the colors
myLabelsColors = feval(fun_definitions);

% Colormap
cmap0  = [0 0 0];
dc  = 0.5;
N = floor(1/dc)+1;
ijk = 1;
for i = 1:N
    for j = 1:N
        for k = 1:N
            cmap(ijk+1,:) = cmap0 + [0 0 k-1] * dc; %#ok<AGROW>
            ijk = ijk + 1;
        end
        cmap0 =  [cmap0(1) j*dc 0];
    end
    cmap0 =  [i*dc 0 0];
end
cnum = 2;

% Colors, labels, line style
lentag = length(TAGS);
mylabels = cell(1,lentag);
mycolors = cell(1,lentag);
mycolors2 = cell(1,lentag);
mylinestyle = cell(1,lentag);
for iT = 1:lentag
    pos = findcellstr(myLabelsColors(:,1),TAGS{iT});
    if pos ~= 0
        mylabels{iT} = myLabelsColors{pos,2};
        mycolors{iT}  = myLabelsColors{pos,3};
        mycolors2{iT}  = myLabelsColors{pos,4};
        if isempty(mycolors2{iT})
            mycolors2{iT}='none';
        end
        mylinestyle{iT}  = myLabelsColors{pos,5};
    else   % defaults
        mylabels{iT}  = TAGS{iT};     % label is just the tag
        mycolors{iT}  = cmap(cnum,:); % graded color map
        mycolors2{iT} = [1 0 0];            % nothing
        mylinestyle{iT}  = '-';          % line
        cnum = cnum + 1;
    end
end   % iT