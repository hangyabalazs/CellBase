function f0=selectivity_plot(time,SEL,params,varargin)
% 


% define_fun = @defineLabelsColors;

default_args = params;
default_args.TitleStr = '';
default_args.XLabel = 'Time';
default_args.YLabel = 'Neurons';
default_args.Thres = [-2 2];
[par, error] = parse_args(default_args,varargin{:});

%----------------------


%---------------------------------
%  Color
%---------------------------------
I1=0.99;
I2=1;
I=I1:(I2-I1)/31:I2;
INT = repmat(I',1,3);
bkgnd=[0 0 0];
map1=make2color([1 0 0],bkgnd,32);
map2=make2color(bkgnd,[0 1 0],32);
cmap = [map1.*INT; map2.*INT(end:-1:1,:)];
%---------------------------------

pos=restrict(time,par.NormalizationWindow(1),par.NormalizationWindow(2));
[D,ind]=sort(sum(SEL(:,pos),2));
ind2=union(find(D<par.Thres(1)),find(D>par.Thres(2)));
ind = ind(ind2);

%
vNULL=[0.5 0];
pNULL=nearest(nanmean(SEL(:)),vNULL);
NULL=vNULL(pNULL);
%
SEL(find(isnan(SEL)))=NULL;
pNEG=find(D(ind2)<NULL);
pPOS=find(D(ind2)>NULL);
%
if   NULL==0
    index = [ind(pPOS); ind(pNEG)];
else
    index = [ind(pPOS(end:-1:1)); ind(pNEG(end:-1:1))];
end
% clear index2
% for i=1:size(S,1)
%      pN(i)=min([find(S(i,:)<thres2(1)) 1000]);
%      pP(i)=min([find(S(i,:)>thres2(2)) 1000]);
%      if pN(i) < pP(i)
%          index2{1,i} = pN(i);
%      else
%          index2{2,i} = pP(i);
%      end
% end
% index2 = [setdiff([index2{1,:}],1000) setdiff([index2{2,:}],1000)];

Dneg=nanmean(SEL(ind(pNEG),:));
Dpos=nanmean(SEL(ind(pPOS),:));
    
%set NaNs to the chance level
NULL=[0.5 0];
pNULL=nearest(nanmean(SEL(:)),NULL);
SEL(find(isnan(SEL)))=NULL(pNULL);

%------------------------------------------------------------
pos2disp=restrict(time,par.window(1),par.window(2));
%alim=[par.window(1) par.window(2) MN MX];
%-------------------------------------------------------------

 
f0=figure(par.FigureNum);
clf; 
imagesc(time(pos2disp),1:length(index),SEL(index,pos2disp));
colormap(cmap);
colorbar;
yl=ylabel(par.YLabel);
xl=xlabel(par.XLabel);
if ~isempty(par.TitleStr)
    tl=title(par.TitleStr);
else
    tl=xl;
end
%axis(alim);
setmyplot(gca,[xl yl tl]);

if strcmpi(par.PlotAverage,'on')
    ax1=gca;
    ax2=axes;
    set(ax2,'Position',get(ax1','Position'));
    hold on;
    if ~isempty(pNEG)
        plot(time(pos2disp),Dneg(pos2disp),'w','LineWidth',2);
    end
    plot(time(pos2disp),Dpos(pos2disp),'w','LineWidth',2);
    set(ax2,'Color','none','YaxisLocation','right');
    set(ax2,'Xticklabel',[])
end
