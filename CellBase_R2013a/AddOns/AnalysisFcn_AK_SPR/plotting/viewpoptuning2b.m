function [NORMstimrateC, NORMstimrateE] = viewpoptuning2b(cellids,varargin) 


if (nargin < 1)
	help viewpoptuning
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
try, g.EpochName;       catch, g.EpochName     = 'AfterChoice'; end;
try, g.Bins;            catch, g.Bins          = 7;   end;
try, g.MinPoints;       catch, g.MinPoints     = 8;  end;
try, g.Plot;            catch, g.Plot          = {'CAF','StimRate'}; end; %{'CAF','CAFdiff','StimRate'}; 
try, g.FigureNum;       catch, g.FigureNum     = 1; end;
try, g.OdorPairID;      catch, g.OdorPairID    = 1; end;
try, g.Pthreshold;      catch, g.Pthreshold    = 0.05; end;
try, g.StimAlign;       catch, g.StimAlign     = 1; end;
try; g.StimNormalization; catch, g.StimNormalization = 'max'; end;  %'mean', 'none', 'sidemean'
try; g.StimCenter;      catch,  g.StimCenter   = 1; end;
try; g.NormCAF;         catch,  g.NormCAF      = 'linear'; end;
    
% test arguments for consistency
% switch lower(g.Plot)
%     case { 'caf', 'caffdiff', 'stimrate'}, ;
%     otherwise error('PLOT must be either on or off');
% end;
FIF = 0;


Xunc=0:0.08:1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumCells = length(cellids);
NORMstimrateC = zeros(NumCells,8)*NaN;
NORMstimrateE = zeros(NumCells,8)*NaN;

fprintf('%d cells to process.\n',NumCells)
for iC = 1:NumCells
 
print_progress(iC,round(NumCells/100),5);

cellid = cellids{iC};
TE = loadcb(cellid,'Events');
SP = loadcb(cellid,'EVENTSPIKES');


epoch_pos = findcellstr(SP.epochs(:,1),g.EpochName);
if (epoch_pos == 0)
  error('Epoch name not found');
end

EpochRate = SP.epoch_rates{epoch_pos};

%valid_trials = find(TE.OdorPokeValid & TE.WaterPokeValid);
%Nmixes = length(unique(TE.OdorRatio(find(TE.OdorPairID == OdorPairID))));
  
if g.OdorPairID == 0    
    trialsC = selecttrial(TE,'OdorConc == 100 & Correct ');
    trialsE = selecttrial(TE,'OdorConc == 100 & Error ');
    trialsL = selecttrial(TE,'OdorConc == 100 & ChoiceLeft');
    trialsR = selecttrial(TE,'OdorConc == 100 & ChoiceRight');
else
    trialsC = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & Correct ',g.OdorPairID));
    trialsE = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & Error ',g.OdorPairID));
    trialsL = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & ChoiceLeft ',g.OdorPairID));
    trialsR = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & ChoiceRight ',g.OdorPairID));
    trialsALL = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100',g.OdorPairID));
end
trials_diff = findpos(TE.OdorRatio,[32 44 56 68]);
trialsCdiff = intersect(trials_diff,trialsC);
trialsEdiff = intersect(trials_diff,trialsE);

if sum(ismember(g.Plot,'StimRate'))
    
[stimC, FstimC, FstimC_se] = tuning_curve(EpochRate(trialsC),TE.OdorRatio(trialsC));
[stimE, FstimE, FstimE_se] = tuning_curve(EpochRate(trialsALL),TE.OdorRatio(trialsALL));

[stimC, FstimC, FstimC_se] = tuning_curve(TE.Correct(trialsALL),TE.OdorRatio(trialsALL));
stimE = 0;
FstimE = 0;
FstimE_se = 0;

fifty=find(stimE==50);
if fifty
    stimE(fifty) = [];
    FstimE(fifty) = [];
    FIF=FIF+1;
end

StandardRatios = [0 20 32 44 56 68 80 100];
stimC = standardize(stimC,StandardRatios);
stimE = standardize(stimE,StandardRatios);

Lc=length(FstimC)/2;        
Le=length(FstimE)/2; 

if g.StimCenter
    pos = findpos(stimC,stimE);
    switch length(stimC)
        case 4,  
             if str2num(cellid(2:3)) == 48
                stimC=[2 3 6 7];   %stimC=[1 3 6 8];
             else
                stimC=[2 3 6 7]; 
             end
        case 6,  stimC=[2 3 4 5 6 7];                      
        case 8,  stimC=[1 2 3 4 5 6 7 8];
    end;
    stimE = stimC(pos);
    %stimE = stimC;
end
%  invert tuning curves to align them to max
if g.StimAlign 
 if FstimC(Lc)>FstimC(Lc+1)
     FstimC = FstimC(end:-1:1);
     FstimE = FstimE(end:-1:1);
 end
end


switch lower(g.StimNormalization)
    case 'none'
        NORMfactor = 1;
    case 'mean'
        NORMfactor = nanmean([FstimC FstimE]);
    case 'max'
        NORMfactor = max([FstimC FstimE]);
    case 'sidemean'
        NORMfactor = [nanmean([FstimC(1:Lc) FstimE(1:Lc)]);  nanmean([FstimC(Lc+1:end) FstimE(Lc+1:end)])];
    case 'maxrate'
        NORMfactor = max(EpochRate(union(trialsC,trialsE)));
    case 'maxratecorr'
        NORMfactor = max(EpochRate(trialsC));
    case 'maxrateerr'
        NORMfactor = max(EpochRate(trialsE));
    case 'perc'   
         NORMfactor = prctile(EpochRate(trialsE),90);
    case 'meanrate'
        NORMfactor = nanmean(EpochRate(union(trialsC,trialsE)))
    otherwise
        NORMfactor = 1;
end

if length(NORMfactor) == 2
    NORMstimrateC(iC,stimC) = [FstimC(1:Lc)/NORMfactor(1) FstimC(Lc+1:end)/NORMfactor(2)];
    NORMstimrateE(iC,stimE) = [FstimE(1:Lc)/NORMfactor(1) FstimE(Lc+1:end)/NORMfactor(2)];
else
    NORMstimrateC(iC,stimC) = FstimC/NORMfactor;
    NORMstimrateE(iC,stimE) = FstimE/NORMfactor;
end

end %ismember StimRate


if sum(ismember(g.Plot,'CAF'))
    
    [cafO , binsO , nt] = histratio(EpochRate(trialsC),EpochRate(trialsE),g.Bins,'fixed_bins','minpoints',g.MinPoints); 
    NORMcafO(iC,:)     = normalize(binsO,cafO,Xunc,g.NormCAF);
    
end % ismember CAF



if sum(ismember(g.Plot,'CAFdiff'))
  [cafOdiff, binsOdiff, nt] = histratio(EpochRate(trialsCdiff),EpochRate(trialsEdiff),g.Bins,'fixed_bins','minpoints',g.MinPoints); 
  NORMcafOdiff(iC,:) = normalize(binsOdiff,cafOdiff,Xunc,g.NormCAF);
end %ismember CAFdiff  

%[cafDIR, binsDIR, nt] = histratio(EpochRate(trialsL),EpochRate(trialsR),g.Bins,'fixed_bins','minpoints',g.MinPoints); 


end % iC

if sum(ismember(g.Plot,'CAF'))
    meanNORMcafO     = nanmean(NORMcafO);
    seNORMcafO       = nanstd(NORMcafO)/sqrt(NumCells-1);
end
if sum(ismember(g.Plot,'CAFdiff'))
    meanNORMcafOdiff = nanmean(NORMcafOdiff);
    seNORMcafOdiff   = nanstd(NORMcafOdiff)/sqrt(NumCells-1);
end
if sum(ismember(g.Plot,'StimRate'))
    
    %Somehow infinities sneak and nanmean doesn't handle it
    NORMstimrateC(find(isinf(NORMstimrateC))) = NaN; 
    NORMstimrateE(find(isinf(NORMstimrateE))) = NaN;
    
    meanNORMstimrateC = nanmean(NORMstimrateC);
    seNORMstimrateC   = nanstd(NORMstimrateC)/sqrt(NumCells-1);
    meanNORMstimrateE = nanmean(NORMstimrateE);
    seNORMstimrateE   = nanstd(NORMstimrateE)/sqrt(NumCells-1);
end
%NFSTIME=NFSTIME-NFSTIMC;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if sum(ismember(g.Plot,'StimRate'))
    
StandardRatios = StandardRatios([1 3 4 5 6 8]);
pC=2:7;
pE=3:6;

[i,j]=find(NORMstimrateE>0);
pj=find(j==2);
if length(pj) < 5
    meanNORMstimrateE(2) = NaN;
end
pj=find(j==7);
if length(pj) < 5
    meanNORMstimrateE(7) = NaN;
end
ONOFF{1} = 'off'; ONOFF{2} = 'on';

figure(g.FigureNum)
clf;
hold on;
eC= errorbar(StandardRatios,meanNORMstimrateC(pC),seNORMstimrateC(pC),'t0-g');
p1=plot(StandardRatios,meanNORMstimrateC(pC),'go');
set(p1,'MarkerFaceColor','g');
eE= errorbar(StandardRatios(2:5),meanNORMstimrateE(pE),seNORMstimrateE(pE),'t0-r');
all=[meanNORMstimrateC meanNORMstimrateE];
axis([-5 105 min(all)*0.9 max(all)*1.1]);
set([eE(2) eC(2)],'LineWidth',3);
set(gca,'XTick',StandardRatios)
xlabel('%A stimulus'); ylabel('Normalized rate');
titstr=sprintf('StimAlign %s; Center %s; Normalize %s',ONOFF{g.StimAlign+1},ONOFF{g.StimCenter+1},g.StimNormalization);
title(titstr);
setmyplot(gca);

g.FigureNum = g.FigureNum+1;
end

if sum(ismember(g.Plot,{'CAF','CAFdiff'}))
    
figure(g.FigureNum)
clf;
hold on;

if sum(ismember(g.Plot,{'CAF'}))
    errorshade(Xunc,meanNORMcafO,seNORMcafO,[0.8 0.2 0],2,[0.8 0.2 0])
end

if sum(ismember(g.Plot,{'CAFdiff'}))
    errorshade(Xunc,meanNORMcafOdiff,seNORMcafOdiff,[0 0.2 0.8],2,[0 0.2 0.8])
end
set(gca,'XTick',0:0.2:1,'YTick',0.5:0.1:1);
setmyplot(gca);
xlabel('Normalized rate'); ylabel('Accuracy');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ncaf = normalize(bins,caf,Xunc,method)

%normalize
if size(bins,2) > 1
    bins = bins/max(bins);
    if strcmp(method,'linear')
        pNaN = find(isnan(caf));
        pNaN = setdiff(pNaN,[1 length(caf)]);
        caf(pNaN) = mean([caf(pNaN-1); caf(pNaN+1)]);
        ncaf = interp1(bins,caf,Xunc,'linear');
    else
        ncaf = interp1(bins,caf,Xunc,'nearest');
    end
else
    ncaf = Xunc*NaN;
end

function stim = standardize(stimC,StandardRatios);

for iS = 1:length(stimC)
     stim(iS) = nearest(StandardRatios,stimC(iS));
end