function fig = plot_behavior(te,varargin)


[P PM N]= get_performance(te);
STIM = sq(PM(1,:,1));

if nargin > 1
    figure(varargin{1}); 
else
    fig=figure;
end
clf; 
sp(1) = subplot(431);

hold on;
for iG = 1:N
    PF = P(iG).psycho;
%     NUM{iG}=[P(iG).psycho(:,5)]';
%     RP{iG}=[P(iG).psycho(:,6)]';
%     STIM{iG}=[P(iG).psycho(:,1)]';
    e=errorbar(PF(:,1),PF(:,2),PF(:,4));
    set(e,'color',gcolor(iG))
end
set(gca,'xtick',STIM)
axis([-5 105 -0.05 1.05]);
ylabel('Fraction choice A');
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
sp(2)=subplot(433);
hold on;
binsO=0:0.05:max(te.OdorSampDur);
binsM=0:0.05:max(te.MovementDur);
nO = hist(te.OdorSampDur,binsO);
nM = hist(te.MovementDur,binsM);
pO=nO/sum(nO); pM=nM/sum(nM);
axO=stairs(binsO,pO,'g');
axM=stairs(binsM+0.01,pM,'c');
set([axO axM],'linewidth',2);
axis([0 max(1.05,max([binsO])*0.95) 0 max([pO pM])*1.1]);
lOM=legend('OSD','MT',1);
set(lOM,'pos',get(lOM,'pos')+[0.08 0.01 0 0]);
set(lOM,'FontSize',8,'Box','off');
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
sp(3)=subplot(437);
NUM=sq(PM(:,:,5))';
bN=bar(STIM,NUM,0.9,'grouped');
for iG = 1:N
    set(bN(iG),'facecolor',gcolor(iG));
end
axis([-3*N 100+3*N 0 max(NUM(:))*1.2]);
ylabel('stim #');
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
sp(4)=subplot(439);
WAIT = te.WaterPokeOut(find(te.RewardMissed))-te.WaterPokeIn(find(te.RewardMissed));
NumWait = length(WAIT);
binsW=0:0.15:max(WAIT);
nW = hist(WAIT,binsW);
pW = nW/sum(nW);
axW=stairs(binsW,pW,'m');
set(axW,'linewidth',2);
axis([0 max(WAIT) 0 max([pW])*1.1]);
%lW=legend('NoRewWait',2);
%set(lW,'FontSize',8,'Box','off');
tW=title(sprintf('%d trials w/o reward',NumWait));
set(tW,'VerticalAlignment','top')
box off
xlabel('Time (s)');
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%

sp(5)=subplot(4,3,10);
RP=sq(PM(:,:,6))';
bRP=bar(STIM,RP,0.9,'grouped');
for iG = 1:N
    set(bRP(iG),'facecolor',gcolor(iG));
end
axis([-3*N 100+3*N 0 1]);
set(gca,'xtick',STIM)
xlabel('A% stimulus');
ylabel('Rew P');



set(sp(1),'pos',[0.1 0.6 0.5 0.33]);
set(sp(2),'pos',[0.7 0.6 0.25 0.33]);
set(sp(3),'pos',[0.1 0.4 0.5 0.13]);
set(sp(4),'pos',[0.7 0.2 0.25 0.33]);
set(sp(5),'pos',[0.1 0.2 0.5 0.13]);


NumConc = length(unique(te.OdorConc));
if NumConc > 1
    text(-4,-1,sprintf('%d concentrations averaged.',NumConc));
end