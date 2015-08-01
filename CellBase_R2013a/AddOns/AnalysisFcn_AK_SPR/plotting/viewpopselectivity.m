function viewpopselectivity(cellids,varargin) 

if (nargin < 1)
	help viewtuning
	return
end

% if  validcellid(cellid,{'list'}) ~= 1
%     fprintf('%s is not valid.',cellid);
%     return
%     
% end
% put arguments into a structure
if ~isempty(varargin)
    try, g = struct(varargin{:}); 
    catch, error('Argument error in the {''param'', value} sequence'); end; 
end;

% assign defaults if no arguments passed
try, g.Type;           catch, g.Type            = 'bar'; end;
try, g.Property;       catch, g.Property        = 'Rate'; end;
try, g.Property1;      catch, g.Property1       = '';   end;
try, g.Property2;      catch, g.Property2       = 'Rate';  end;
try, g.Significance;   catch, g.Significance    = []; end;
try, g.Significance1;  catch, g.Significance1   = []; end;
try, g.Significance2;  catch, g.Significance2   = []; end;
try, g.Pthreshold;     catch, g.Pthreshold      = 0.05; end;
try, g.ValueThreshold; catch, g.ValueThreshold  = []; end;
try, g.NumBins;        catch, g.NumBins         = 20; end;
try, g.FigureNum;      catch, g.FigureNum       = 1; end;
try, g.LeftColor;      catch, g.LeftColor       = 'blue'; end;
try, g.RightColor;     catch, g.RightColor      = 'blue'; end;
try  g.ColorBoth;      catch, g.ColorBoth       = [0.5 0 0.5]; end;
try  g.ColorNeither;   catch, g.ColorNeither    = [0.8 0.8 0.8]; end;
try  g.Color1;         catch, g.Color1          = [0.8 0 0.2]; end;
try  g.Color2;         catch, g.Color2          = [0 0.2 0.8]; end;
try  g.Marker;         catch, g.Marker          = 'o'; end;
try  g.Transform1;     catch, g.Transform1      = ''; end;
try  g.Transform2;     catch, g.Transform2      = ''; end;
% test arguments for consistency


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(g.Type,'bar')
    
D=getvalue(g.Property,cellids);
[D, ind]=setdiff(D,[Inf -Inf]);


if ~isempty(g.Significance)
    P=getvalue(g.Significance,cellids);
    P = P(ind);
    Psig=find(P<g.Pthreshold);
else
    g.Pthreshold = 0;
    Psig = 1:length(D);
end

if ~isempty(g.ValueThreshold)
    Dsmall = intersect( find(D<g.ValueThreshold) , find(D>-g.ValueThreshold));
else
    Dsmall = [];
end
GOOD = setdiff(Psig,Dsmall);
NOTSIG = setdiff(1:length(D),GOOD);

dDATA = (max(D)-min(D))/g.NumBins;
edges = min(D):dDATA:max(D)+2*eps;
bins = edges+dDATA;
nG  =histc(D(GOOD),edges);
if ~isempty(NOTSIG)
    nNS =histc(D(NOTSIG),edges);
else
    nNS = nG*0;
end

nG = nG(:);
nNS = nNS(:);
% 
% NEG = GOOD(find(D(GOOD)<0));
% POS = GOOD(find(D(GOOD)>0));
%stat_text = sprintf('%d/%d (p<%.2f)  +%d:%.2f -%d:%.2f',length(GOOD),length(D),g.Pthreshold,length(POS),nanmedian(D(POS)),length(NEG),nanmedian(D(NEG)));

stat_text = sprintf('%d/%d (p<%.2f) median = %.2f',length(GOOD),length(D),g.Pthreshold,nanmedian(D));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure(g.FigureNum)
clf;

b=bar(bins,[nG nNS],'stacked');
set(b(2),'facecolor',[0.8 0.8 0.8]);
axis([bins(1)-dDATA bins(end)+dDATA 0 max([nNS; nG])*1.1])
xlabel(g.Property);
ylabel('Count')
title(stat_text);
setmyplot(gca);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif strcmpi(g.Type,'scatter')
    
D1=getvalue(g.Property1,cellids);
D2=getvalue(g.Property2,cellids);

switch lower(g.Transform1)
    case {'absolute','abs'}
        D1 = abs(D1);
    case {'invert','inv'}
        D1 = 0.5+D1;
end
switch lower(g.Transform2)
    case {'absolute','abs'}
        D2 = abs(D2);
    case {'invert','inv'}
        D2 = 0.5+D2;
end

if ~isempty(g.Significance1)
    P1=getvalue(g.Significance1,cellids);
    P1sig=find(P1<g.Pthreshold);
else
    P1sig = 1:length(D1);
end

if ~isempty(g.Significance2)
    P2=getvalue(g.Significance2,cellids);
    P2sig=find(P2<g.Pthreshold);
else
    P2sig = 1:length(D2);
end


if ~isempty(g.ValueThreshold)
    D1small = intersect( find(D1<g.ValueThreshold) , find(D1>-g.ValueThreshold));
    D2small = intersect( find(D2<g.ValueThreshold) , find(D2>-g.ValueThreshold));
else
    D1small = [];
    D2small = [];
end

D1good = setdiff(P1sig,D1small);
D2good = setdiff(P2sig,D2small);

BOTH = intersect(D1good,D2good);
ONLY1 = setdiff(D1good,BOTH);
ONLY2 = setdiff(D2good,BOTH);
NEITHER = setdiff(1:length(D1),union(D1good,D2good));

[r,p]=corrcoef(D1(BOTH),D2(BOTH));
cc=r(2,1);
ccSIG=(p(2,1)<g.Pthreshold)+1;
STAR = ' *';

stat_text = sprintf('BOTH/%s/%s/NS: %d/%d/%d/%d (p<%.2f)  CC:%.2f%s',g.Property1,g.Property2,length(BOTH),length(ONLY1),length(ONLY2),length(NEITHER),...
            g.Pthreshold,cc,STAR(ccSIG));


figure(g.FigureNum)
clf; hold on
p4=plot(D1(NEITHER),D2(NEITHER),'ko');
set(p4,'Marker',g.Marker,'Color','w','MarkerFaceColor',g.ColorNeither,'MarkerEdgeColor','none');

p2=plot(D1(ONLY1),D2(ONLY1),'ko');
set(p2,'Marker',g.Marker,'Color','w','MarkerFaceColor',g.Color1,'MarkerEdgeColor','none');
p3=plot(D1(ONLY2),D2(ONLY2),'ko');
set(p3,'Marker',g.Marker,'Color','w','MarkerFaceColor',g.Color2,'MarkerEdgeColor','none');
p1=plot(D1(BOTH),D2(BOTH),'ko');
set(p1,'Marker',g.Marker,'Color','w','MarkerFaceColor',g.ColorBoth,'MarkerEdgeColor','none');
l=legend('Neither',g.Property1,g.Property2,'Both');
set(l,'box','off');
xlabel(g.Property1); ylabel(g.Property2);
title(stat_text);
ax=axis;
%axis([-0.05 ax(2) -0.05 ax(4)]);
setmyplot(gca);

end