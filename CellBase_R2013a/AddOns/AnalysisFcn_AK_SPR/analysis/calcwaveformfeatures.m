function out=calcwaveformfeatures(waveforms,varargin)
% Calculates features for all waveforms: new features need to be defined
% INPUTS:
% waveforms: ideally should be able to handle single waveform, multiple
% waveforms and multiple channels... i.e. vector, 2D matrix, 3D matrix
% varargin: ideally I should be able to select which properties to
% calculate. Currently, all will be calculated.
% OUTPUT: structure with properties as fields. Each waveform will have its
% own properties.

% Deal with default arguments
default_args={...
    'props',           'all';... % if you want to retrieve all props, otherwise, mention the props
    'sample_rate',      32*10^3;...
    'peak_index',       7:11;...
    };
[g, error] = parse_args(default_args,varargin{:});

%% DEAL WITH WAVEFORMS TO DETERMINE NATURE OF INPUT

%% waveforms is a (numspikes X chan X datapoints)
Ch1= squeeze(waveforms(:,1,:));
Ch2= squeeze(waveforms(:,2,:));
Ch3= squeeze(waveforms(:,3,:));
Ch4= squeeze(waveforms(:,4,:));

clear waveform
% Find largest channel
meanWV =[nanmean(Ch1); nanmean(Ch2); nanmean(Ch3); nanmean(Ch4)];
meanWV;
[YY,II]=max(meanWV(:,8));% this may not be the largest channel since

switch II
    case 1
        n = size(Ch1(:,1));
    case 2
        n = size(Ch2(:,1));
    case 3
        n = size(Ch3(:,1));
    case 4
        n = size(Ch4(:,1));
end

% Align spikes to peak
%         eval(['n=size(Ch' num2str(II) '(:,1));'])
WV=repmat(NaN,n,96);

for i=1:n
    eval(['[YYY,I]=max(Ch' num2str(II) '(i,:));'])
    shift=40-I;
    eval(['WV(i,shift:shift+31)=Ch' num2str(II) '(i,:);'])
end
clear Ch1 Ch2 Ch3 Ch4
% Smooth Waveform
X =0:32;
XX =0:.1:32;
Y = squeeze(nanmedian(WV));
SY=squeeze(nanstd(WV));
yy =spline(X,Y(31:63),XX);
syy =spline(X,SY(31:63),XX);
% clear WV
% HACK
if isnan(mean(yy))
    yy = spline(X,Y(32:64),XX);
end
% Find Peak Value and Position
[p,pi] =max(yy);
Peak = p;
PeakI= XX(pi);
% PeakI= pi;

% Find Post-Vally value and position
[v,vi]= min(yy(pi:end));
Valley = v;
ValleyI =XX(vi+pi);
% ValleyI =vi+pi;

% Find Pre-Valley value and position
[pv,pvi] =min(yy(1:pi));
PreValley = pv;
PreValleyI =XX(pvi);

% Find width of hyperpolarization
base = median(yy(1:40)); % arbitrarily designated as baseline
wsi = find(yy(pi:vi+pi)<=base, 1 );
wsi = (pi-1)+wsi;
wei = (vi+pi-1)+find(yy(vi+pi:end)>=base, 1 );
HyperBeginI = XX(wsi);
if ~isempty(wei),HyperEndI = XX(wei);else HyperEndI = XX(end);end % in case it does not return to baseline

dt = 1/g.sample_rate * 1000 / 10; % (per/sec)*100

PV2 = abs(Peak/Valley);
V1V2 =abs(PreValley/Valley);
Width =(ValleyI-PeakI)*dt*10000; % (per/sec)*100*10000 i.e. microsec
Hyper = (HyperEndI-HyperBeginI)*dt*10000;
Width2 = (HyperEndI-PeakI)*dt*10000;
Repolar =(Peak-Valley)/Width;
% y=1;
% Data(y,1)={'Cell_path'};
% Data(y,2)={cellid};y=y+1;
% Data(y,1)={'WaveformI'};
% Data(y,2)={XX};y=y+1;
% Data(y,1)={'Waveform'};
% Data(y,2)={yy};y=y+1;
% Data(y,1)={'Std Waveform'};
% Data(y,2)={syy};y=y+1;
% Data(y,1)={'Peak'};
% Data(y,2)={Peak};y=y+1;
% Data(y,1)={'Peak_location'};
% Data(y,2)={PeakI};y=y+1;
% Data(y,1)={'Valley1'};
% Data(y,2)={PreValley};y=y+1;
% Data(y,1)={'Valley1_location'};
% Data(y,2)={PreValleyI};y=y+1;
% Data(y,1)={'Valley2'};
% Data(y,2)={Valley};y=y+1;
% Data(y,1)={'Valley2_location'};
% Data(y,2)={ValleyI};y=y+1;
% Data(y,1)={'Width'};
% Data(y,2)={Width};y=y+1;
% Data(y,1)={'PV2_Ratio'};
% Data(y,2)={PV2};y=y+1;
% Data(y,1)={'Valleys_ratio'};
% Data(y,2)={V1V2};y=y+1;
% Data(y,1)={'Repolarization'};
% Data(y,2)={Repolar};y=y+1;
% Data(y,1)={'Hyperpolarization Index'};
% Data(y,2)={[HyperBeginI HyperEndI]};y=y+1;
% Data(y,1)={'Hyperpolarization'};
% Data(y,2)={Hyper};y=y+1;
% Data(y,1)={'Width2'};
% Data(y,2)={(HyperEndI-PeakI)*dt*10000};y=y+1;
% Data(y,1)={'Spike Times sec'};
% Data(y,2)={SpikeTimes*1e-4};y=y+1;

Data=struct;
Data.cellid=cellid;
Data.WaveformI=XX;
Data.Waveform=yy;
Data.WaveformStd=syy;
Data.Peak=Peak;
Data.Peak_location=PeakI;
Data.Valley1=PreValley;
Data.Valley1_location=PreValleyI;
Data.Valley2=Valley;
Data.Valley2_location=ValleyI;
Data.Width=Width;
Data.PV2_Ratio=PV2;
Data.Valleys_ratio=V1V2;
Data.Repolarization=Repolar;
Data.HyperpolarizationI=[HyperBeginI HyperEndI];
Data.Hyperpolarization=Hyper;
Data.Width2=(HyperEndI-PeakI)*dt*10000;
Data.stimes=SpikeTimes;
Data.waveforms=wave;


