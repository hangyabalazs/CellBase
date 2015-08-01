% CALC_LIGHT_SPONTANEOUS_WAVEFORM_FEATURES
% selects spontaneous and light activated spikes, extracts the respective
% waveforms and calculates waveform features and stores it in
% SPONT_WAVEDATA AND LIGHT_WAVEDATA that are saved on the cellbase datapath
% along with other parameters used.
% SPR 2010-10-29

SPONT_WAVEDATA=struct;
LIGHT_WAVEDATA=struct;
tic;
problem_waveform_cellids=[];
fields={'cellid','WaveformI','Waveform','WaveformStd','Peak',...
    'Valley2','Amplitude','Peak_location','Peak_fwhh','Peak_fwqh',...
    'Valley_fwhh','Valley1','Valley1_location','Valley2_location',...
    'Width','RestVoltage','PV2_Ratio','Valleys_ratio','Repolarization',...
    'HyperpolarizationI','Hyperpolarization','Width2','stimes','waveforms',...
    'meanrate'};
% INITIALIZATION
for iF=1:length(fields),
    SPONT_WAVEDATA(1).(fields{iF})=[];
end
LIGHT_WAVEDATA=SPONT_WAVEDATA;
% cellids=cellids(876:end);
cellids=listtag('cells');
wfparams.margins=[0.00 0.004];
wfparams.MaxSpikes=10000;
wfparams.whichpulses='BurstOn'; % PulseOn
wfparams.ifsave=0;
lw=2; % linewidth for online plotting
for iCell=1:length(cellids),
    cellid=cellids(iCell);
    try
        EV=loadcb(cellid,'Events');
        try    SE=loadcb(cellid,'StimEvent'); catch, SE=[]; end
        try    TE=loadcb(cellid,'TrialEvent'); catch, TE=[]; end
        
        % exclude light stimulation epochs
        try
            %     fp=find(SE.FirstPulse==1);
            bon=SE.BurstOn(~isnan(SE.BurstOn));
            boff=SE.BurstOff(~isnan(SE.BurstOff));
            sess_start=EV.Events_TimeStamps(1);
            sess_end=EV.Events_TimeStamps(end);
            stpulses=[bon boff sess_start sess_end];
        catch
            stpulses=[];
        end
        
        % exclude light stimulation epochs during behavior
        try
            evpulses=[TE.AllPulseOns{~isempty(TE.AllPulseOns)}+wfparams.margins(1) TE.AllPulseOns{~isempty(TE.AllPulseOns)}+wfparams.margins(2)...
                TE.PulseOn(~isnan(TE.PulseOn))+wfparams.margins(1) TE.PulseOn(~isnan(TE.PulseOn))+wfparams.margins(2)];
        catch
            evpulses=[];
        end
        allpulses=sort([stpulses evpulses]);
        tseg=reshape(allpulses,2,[]);

        % if there are any time segments with spontaneous activity
        if ~isempty(tseg),
            [nonstimts,nonstimind,nonstimisi]=extractSegSpikes(cellid,tseg);
            % if we have too many spikes lets take the MaxNumSpikes
            if length(nonstimts)>wfparams.MaxSpikes,
                rand_i=randperm(length(nonstimts));
                nonstimts=nonstimts(rand_i(1:wfparams.MaxSpikes));
            end
            % if we got any spontaneous spikes, lets calculate spontaneous
            % waveform features
            if length(nonstimts)>1,
                [spont_wavefeatures,junk]=get_mclust_waveforms4(cellid,nonstimts*1e4,[]);
                spont_wavefeatures.meanrate=1/nanmean(nonstimisi);
                SPONT_WAVEDATA(iCell)=spont_wavefeatures;
            elseif length(nonstimts)==1,
                [spont_wavefeatures,junk]=get_mclust_waveforms4(cellid,nonstimts*1e4,[]);
                spont_wavefeatures.meanrate=NaN;
                SPONT_WAVEDATA(iCell)=spont_wavefeatures;
            else
                SPONT_WAVEDATA(iCell).cellid=cellid;
            end
        else % no time segments, no spikes
            SPONT_WAVEDATA(iCell).cellid=cellid;
        end
        
        % Now lets get the light activated spikes
        try
            switch wfparams.whichpulses
                case 'BurstOn'
                    fp=find(SE.FirstPulse==1);
                case 'PulseOn'
                    fp=find(~isnan(SE.PulseOn));
            end
            stpulses=SE.PulseOn(fp);
            stsegs=[stpulses+wfparams.margins(1);stpulses+wfparams.margins(2)];
        catch
            stsegs=[];
        end
        try
            evpulses=[TE.AllPulseOns{~isempty(TE.AllPulseOns)} TE.PulseOn(~isnan(TE.PulseOn))];
            evsegs=[evpulses+wfparams.margins(1);evpulses+wfparams.margins(2)];
        catch
            evsegs=[];
        end
        tseg=[stsegs evsegs];
        
        % if there are any periods of light activation
        if ~isempty(tseg),
            [stimts,stimind,stimisi]=extractSegSpikes(cellid,tseg);
            % if we got any light evoked spikes, lets get the light activated waveforms and
            % their features
            if length(stimts)>1,
                [light_wavefeatures,junk]=get_mclust_waveforms4(cellid,stimts*1e4,[]);
                light_wavefeatures.meanrate=1/nanmean(stimisi);
                LIGHT_WAVEDATA(iCell)=light_wavefeatures;
            elseif length(stimts)==1,
                [light_wavefeatures,junk]=get_mclust_waveforms4(cellid,stimts*1e4,[]);
                light_wavefeatures.meanrate=NaN;
                LIGHT_WAVEDATA(iCell)=light_wavefeatures;
            else
                LIGHT_WAVEDATA(iCell).cellid=cellid;
            end
        else % no light activated pertiods, no spikes
            LIGHT_WAVEDATA(iCell).cellid=cellid;
        end
    catch
        problem_waveform_cellids=[problem_waveform_cellids cellid];
    end
    figure(103)
    clf
    subplot(211)
    if size(SPONT_WAVEDATA(iCell).waveforms,1)>1,
        try plot(linspace(0,1,size(SPONT_WAVEDATA(iCell).waveforms,2)),nanmean(SPONT_WAVEDATA(iCell).waveforms,1),'g','LineWidth',lw); catch end
    else
        try plot(linspace(0,1,size(SPONT_WAVEDATA(iCell).waveforms,2)),SPONT_WAVEDATA(iCell).waveforms,'g','LineWidth',lw); catch end
    end
    hold on
    if size(LIGHT_WAVEDATA(iCell).waveforms,2)>1,
        try plot(linspace(0,1,size(LIGHT_WAVEDATA(iCell).waveforms,2)),nanmean(LIGHT_WAVEDATA(iCell).waveforms,1),'c','LineWidth',lw); catch end
    else
        try plot(linspace(0,1,size(LIGHT_WAVEDATA(iCell).waveforms,1)),LIGHT_WAVEDATA(iCell).waveforms,'c','LineWidth',lw); catch end
    end
    
    subplot(212)
    try plot(SPONT_WAVEDATA(iCell).WaveformI,SPONT_WAVEDATA(iCell).Waveform,'g','LineWidth',lw); catch end
    hold on
    try plot(LIGHT_WAVEDATA(iCell).WaveformI,LIGHT_WAVEDATA(iCell).Waveform,'c','LineWidth',lw); catch end
    
    try line([SPONT_WAVEDATA(iCell).Peak_location SPONT_WAVEDATA(iCell).Peak_location],ylim,'Color','k'); catch end
    try line([SPONT_WAVEDATA(iCell).Valley2_location SPONT_WAVEDATA(iCell).Valley2_location],ylim,'Color','k'); catch end
    
    fstamp(cellid,10,'top-center')
end
analysis_timestamp=datestr(now);
if wfparams.ifsave==1,
    save([getpref('cellbase','datapath') filesep 'WAVEDATA'],'wfparams','cellids','LIGHT_WAVEDATA','SPONT_WAVEDATA','analysis_timestamp')
end
toc;