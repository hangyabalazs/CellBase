function TE = solo2trialevents_auditory_pavlovian2(filepath,ifsave)
%SOLO2TRIALEVENTS_AUDITORY_PAVLOVIAN2   Create trial events structure.
%   SOLO2TRIALEVENTS_AUDITORY_PAVLOVIAN2(SESSPATH) creates trial events
%   structure ('TrialEvents') for the session given in SESSPATH. The event
%   are extracted from the saved behavior file. The first trial is omitted.
%   Most of the events are given relative to the TrialStart event.
%
%   SOLO2TRIALEVENTS_AUDITORY_PAVLOVIAN2(SESSPATH,IFSAVE) uses an addidional
%   logical input argument do determine whether to save the results.
%
%   See also MAKETRIALEVENTS2_GONOGO and DEFINEEVENTSEPOCHS_PAVLOVIAN.

%   Edit log: BH 5/13/14

% We want to save the result
if nargin < 2
    ifsave = 1;
end

% Determine filepath and load data
warning off
load(filepath)
warning backtrace

% Session ID (date and 'a' or 'b') implement solofilename2tags here
sessionID = filepath(length(filepath)-10:length(filepath)-4);
sess_datetime = saved.SavingSection_SaveTime;

% No. of trials (exclude trial 1 and last trial may not have been completed
ntrials = saved.ProtocolsSection_n_completed_trials - 1;

% Preallocation of space
[TE.R1Trial TE.R2Trial TE.RTrial TE.P1Trial TE.P2Trial TE.PTrial ...
    TE.O1Trial TE.O2Trial TE.OTrial] = deal(nan(1,ntrials));
[TE.TrialStart TE.TrialEnd TE.StimulusOn TE.StimulusOff TE.DelayStart TE.DelayEnd ...
    TE.Reward1On TE.Reward1Off TE.Reward2On TE.Reward2Off TE.RewardOn TE.RewardOff ...
    TE.Punish1On TE.Punish1Off TE.Punish2On TE.Punish2Off TE.PunishOn TE.PunishOff ...
    TE.Reward1 TE.Reward2 TE.Reward TE.Punish1 TE.Punish2 TE.Punish ...
    TE.DeliverFeedback TE.DeliverAllFeedback ...
    TE.Omission1 TE.Omission2 TE.Omission ...
    TE.RestTimeStart TE.RestTimeEnd] = deal(nan(1,ntrials));
[TE.DelayLick TE.FeedbackLick TE.FirstLick] = deal(nan(1,ntrials));
[TE.ReactionTime TE.Reward1RT TE.Reward2RT TE.RewardRT ...
    TE.Punish1RT TE.Punish2RT TE.PunishRT ...
    TE.Omission1RT TE.Omission2RT] = deal(nan(1,ntrials));
[TE.ITIDistribution TE.TotalITI TE.FirstITIBegins] = deal(nan(1,ntrials));
[TE.LickIn TE.LickOut TE.sessionID TE.datetime ...
    TE.ITIBegins TE.ITIEnds] = deal(cell(1,ntrials));
[TE.LeftWaterValveDur TE.LeftWaterValveOn TE.LeftWaterVolume] = deal(nan(1,ntrials));

% Correct side
TE.RewardBiasTrials = saved.auditory_pavlovian2_TrialTypeList(1:ntrials) == 1;
TE.RewardBiasTrials = [TE.RewardBiasTrials(2:end) NaN];  % remove first trial
TE.PunishBiasTrials = saved.auditory_pavlovian2_TrialTypeList(1:ntrials) == 0;
TE.PunishBiasTrials = [TE.PunishBiasTrials(2:end) NaN];  % remove first trial

% Type of trial (for e.g. diff stim durations get different index)
TE.StimulusName = saved.auditory_pavlovian2_SoundStim(1:ntrials);
TE.StimulusName = [TE.StimulusName(2:end) NaN];  % remove first trial
TE.StimulusID = nan(1,ntrials);
TE.StimulusID(strcmp(TE.StimulusName,'RewardBias')) = 1;   % high prob. reward, low prob. punishment
TE.StimulusID(strcmp(TE.StimulusName,'PunishBias')) = 2;   % high prob. punishment, low prob. reward

TE.FeedbackName = saved.auditory_pavlovian2_Feedback(1:ntrials);
TE.FeedbackName = [TE.FeedbackName(2:end) NaN];  % remove first trial
TE.FeedbackID = nan(1,ntrials);
TE.FeedbackID(strcmp(TE.FeedbackName,'Reward1')) = 1;   % high prob. reward (rewardbias)
TE.FeedbackID(strcmp(TE.FeedbackName,'Reward2')) = 2;   % low prob. reward (punishbias)
TE.FeedbackID(strcmp(TE.FeedbackName,'Punish1')) = 3;   % low prob. punishment (rewardbias)
TE.FeedbackID(strcmp(TE.FeedbackName,'Punish2')) = 4;   % high prob. punishment (punishbias)
TE.FeedbackID(strcmp(TE.FeedbackName,'Omission')) = 5;   % omission

% Animal name
TE.Ratname = saved_history.SavingSection_ratname(1:ntrials)';
TE.Ratname = [TE.Ratname(2:end) NaN];  % remove first trial

% Define variables trialwise
countr = 0;
for iT = 2:ntrials
    countr = countr + 1;
    
    % Trial type
    TE.R1Trial(countr) = zero2nan(double(isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'reward1')));
    TE.R2Trial(countr) = zero2nan(double(isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'reward2')));
    TE.RTrial(countr) = zero2nan(isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'reward1')+...
        isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'reward2'));
    TE.P1Trial(countr) = zero2nan(double(isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'punish1')));
    TE.P2Trial(countr) = zero2nan(double(isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'punish2')));
    TE.PTrial(countr) = zero2nan(isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'punish1')+...
        isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'punish2'));
    TE.O1Trial(countr) = zero2nan(isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'omission')*...
        strcmp(TE.StimulusName{countr},'RewardBias'));   % high prob. reward, omission
    TE.O2Trial(countr) = zero2nan(isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'omission')*...
        strcmp(TE.StimulusName{countr},'PunishBias'));   % low prob. reward, omission
    TE.OTrial(countr) = zero2nan(double(isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'omission')));   % all omissions
    
    % State time stamps
    TE.TrialStart(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.state_0(1,2);
    TE.TrialEnd(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.state_0(2,1);
    
    if ~isempty(saved_history.ProtocolsSection_parsed_events{iT}.states.stimulus)   % stimulus state is empty for Operant trials
        TE.StimulusOn(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.stimulus(1);
        TE.StimulusOff(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.stimulus(2);
        switch TE.FeedbackName{countr}
            case 'Reward1'   % high prob. reward
                TE.Reward1On(countr) = TE.StimulusOn(countr);
                TE.Reward1Off(countr) = TE.StimulusOff(countr);
                TE.RewardOn(countr) = TE.StimulusOn(countr);
                TE.RewardOff(countr) = TE.StimulusOff(countr);
            case 'Reward2'   % low prob. reward
                TE.Reward2On(countr) = TE.StimulusOn(countr);
                TE.Reward2Off(countr) = TE.StimulusOff(countr);
                TE.RewardOn(countr) = TE.StimulusOn(countr);
                TE.RewardOff(countr) = TE.StimulusOff(countr);
            case 'Punish1'   % low prob. punishment
                TE.Punish1On(countr) = TE.StimulusOn(countr);
                TE.Punish1Off(countr) = TE.StimulusOff(countr);
                TE.PunishOn(countr) = TE.StimulusOn(countr);
                TE.PunishOff(countr) = TE.StimulusOff(countr);
            case 'Punish2'   % high prob. punishment
                TE.Punish2On(countr) = TE.StimulusOn(countr);
                TE.Punish2Off(countr) = TE.StimulusOff(countr);
                TE.PunishOn(countr) = TE.StimulusOn(countr);
                TE.PunishOff(countr) = TE.StimulusOff(countr);
        end
        
        TE.DelayStart(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.delay(1);
        TE.DelayEnd(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.delay(2);
        
        if isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'reward1') && ...
                ~isempty(saved_history.ProtocolsSection_parsed_events{iT}.states.reward1)
            TE.Reward1(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.reward1(1);
            TE.Reward(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.reward1(1);
            TE.DeliverFeedback(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.reward1(1);
            TE.DeliverAllFeedback(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.reward1(1);
        end
        if isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'reward2') && ...
                ~isempty(saved_history.ProtocolsSection_parsed_events{iT}.states.reward2)
            TE.Reward2(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.reward2(1);
            TE.Reward(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.reward2(1);
            TE.DeliverFeedback(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.reward2(1);
            TE.DeliverAllFeedback(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.reward2(1);
        end
        if isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'punish1') && ...
                ~isempty(saved_history.ProtocolsSection_parsed_events{iT}.states.punish1)
            TE.Punish1(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.punish1(1);
            TE.Punish(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.punish1(1);
            TE.DeliverFeedback(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.punish1(1);
            TE.DeliverAllFeedback(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.punish1(1);
        end
        if isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'punish2') && ...
                ~isempty(saved_history.ProtocolsSection_parsed_events{iT}.states.punish2)
            TE.Punish2(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.punish2(1);
            TE.Punish(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.punish2(1);
            TE.DeliverFeedback(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.punish2(1);
            TE.DeliverAllFeedback(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.punish2(1);
        end
        if isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'omission')
            switch TE.StimulusName{countr}
                case 'RewardBias'   % high prob. reward
                    TE.Omission1(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.omission(1);
                case 'PunishBias'   % high prob. punishment
                    TE.Omission2(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.omission(1);
            end
            TE.Omission(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.omission(1);
            TE.DeliverAllFeedback(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.omission(1);
        end
        
        TE.RestTimeStart(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.resttime(1);
        TE.RestTimeEnd(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.resttime(2);
    end
    
    % Response variables
    TE.LickIn{countr} = saved_history.ProtocolsSection_parsed_events{iT}.pokes.L(:,1);
    TE.LickOut{countr} = saved_history.ProtocolsSection_parsed_events{iT}.pokes.L(:,2);
    licks = TE.LickIn{countr};
    if ~isempty(licks(licks>TE.DelayStart(countr)&licks<TE.DelayEnd(countr)))
        TE.DelayLick(countr) = 1;
    end
    if ~isempty(licks(licks>TE.DeliverAllFeedback(countr)&licks<TE.RestTimeEnd(countr)))
        TE.FeedbackLick(countr) = 1;
        TE.FirstLick(countr) = licks(find(licks>TE.DeliverAllFeedback(countr),1,'first'));
    end
    
    TE.ReactionTime(countr) = TE.FirstLick(countr) - TE.DeliverAllFeedback(countr);
    if ~isfield(saved_history.ProtocolsSection_parsed_events{iT}.states,'omission')
        switch TE.FeedbackName{countr}
            case 'Reward1'
                TE.Reward1RT(countr) = TE.ReactionTime(countr);
                TE.RewardRT(countr) = TE.ReactionTime(countr);
            case 'Reward2'
                TE.Reward2RT(countr) = TE.ReactionTime(countr);
                TE.RewardRT(countr) = TE.ReactionTime(countr);
            case 'Punish1'
                TE.Punish1RT(countr) = TE.ReactionTime(countr);
                TE.PunishRT(countr) = TE.ReactionTime(countr);
            case 'Punish2'
                TE.Punish2RT(countr) = TE.ReactionTime(countr);
                TE.PunishRT(countr) = TE.ReactionTime(countr);
        end
    else
        switch TE.StimulusName{countr}
            case 'RewardBias'
                TE.Omission1RT(countr) = TE.ReactionTime(countr);
            case 'PunishBias'
                TE.Omission2RT(countr) = TE.ReactionTime(countr);
        end
    end
    
    % ITI variables
    TE.ITIDistribution(countr) = saved.auditory_pavlovian2_ITIs(iT);
    TE.ITIBegins{countr} = saved_history.ProtocolsSection_parsed_events{iT}.states.iti(:,1);
    TE.ITIEnds{countr} = saved_history.ProtocolsSection_parsed_events{iT}.states.iti(:,2);
    TE.FirstITIBegins(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.iti(1,1);
    TE.TotalITI(countr) = TE.ITIEnds{countr}(end) - TE.ITIBegins{countr}(1);
    
    % Valve variables
    TE.LeftWaterValveOn(countr) = TE.Reward(countr);
    TE.LeftWaterVolume(countr) = saved_history.WaterValvesSection_Left_volume{iT};
    TE.LeftWaterValveDur(countr) = saved.auditory_pavlovian2_ValveOpenTime(iT);
    
    % Session ID
    TE.sessionID{countr} = sessionID;
    TE.datetime{countr} = sess_datetime;
end

% Sound duration
TE.SoundDuration = nan(1,ntrials);
TE.SoundDuration(:) = saved.auditory_pavlovian2_SoundDuration;
TE.SoundDuration = [TE.SoundDuration(2:end) NaN];  % remove first trial

% Sound Frequency
TE.SoundFrequency = cell(1,ntrials);
TE.SoundFrequency(:) = {saved.auditory_pavlovian2_SoundFrequency};
TE.SoundFrequency = [TE.SoundFrequency(2:end) NaN];  % remove first trial

% Convert to relative time to TrialStart
TE.TrialEnd = TE.TrialEnd - TE.TrialStart;
TE.StimulusOn = TE.StimulusOn - TE.TrialStart;
TE.StimulusOff = TE.StimulusOff - TE.TrialStart;
TE.Reward1On = TE.Reward1On - TE.TrialStart;
TE.Reward1Off = TE.Reward1Off - TE.TrialStart;
TE.Reward2On = TE.Reward2On - TE.TrialStart;
TE.Reward2Off = TE.Reward2Off - TE.TrialStart;
TE.Punish1On = TE.Punish1On - TE.TrialStart;
TE.Punish1Off = TE.Punish1Off - TE.TrialStart;
TE.Punish2On = TE.Punish2On - TE.TrialStart;
TE.Punish2Off = TE.Punish2Off - TE.TrialStart;
TE.RewardOn = TE.RewardOn - TE.TrialStart;
TE.RewardOff = TE.RewardOff - TE.TrialStart;
TE.PunishOn = TE.PunishOn - TE.TrialStart;
TE.PunishOff = TE.PunishOff - TE.TrialStart;
TE.DelayStart = TE.DelayStart - TE.TrialStart;
TE.DelayEnd = TE.DelayEnd - TE.TrialStart;
TE.Reward1 = TE.Reward1 - TE.TrialStart;
TE.Reward2 = TE.Reward2 - TE.TrialStart;
TE.Reward = TE.Reward - TE.TrialStart;
TE.Punish1 = TE.Punish1 - TE.TrialStart;
TE.Punish2 = TE.Punish2 - TE.TrialStart;
TE.Punish = TE.Punish - TE.TrialStart;
TE.DeliverFeedback = TE.DeliverFeedback - TE.TrialStart;
TE.DeliverAllFeedback = TE.DeliverAllFeedback - TE.TrialStart;
TE.Omission1 = TE.Omission1 - TE.TrialStart;
TE.Omission2 = TE.Omission2 - TE.TrialStart;
TE.Omission = TE.Omission - TE.TrialStart;
TE.RestTimeStart = TE.RestTimeStart - TE.TrialStart;
TE.RestTimeEnd = TE.RestTimeEnd - TE.TrialStart;
TE.FirstLick = TE.FirstLick - TE.TrialStart;
TE.LeftWaterValveOn = TE.LeftWaterValveOn - TE.TrialStart;
TE.FirstITIBegins = TE.FirstITIBegins - TE.TrialStart;
for k = 1:ntrials
    TE.ITIBegins{k} = TE.ITIBegins{k} - TE.TrialStart(k);   % make timestamps relative to TrialStart
    TE.ITIEnds{k} = TE.ITIEnds{k} - TE.TrialStart(k);
end

wasnan = false;   % licks
for k = 1:ntrials
    TE.LickIn{k} = TE.LickIn{k} - TE.TrialStart(k);   % make timestamps relative to TrialStart
    TE.LickOut{k} = TE.LickOut{k} - TE.TrialStart(k);
    naninxin = isnan(TE.LickIn{k});    % handle NaNs (assuming they correspond to trial switches with the animal poked in)
    naninxout = isnan(TE.LickOut{k});
    if any(naninxin) || any(naninxout)
        if ~isequal(sum(naninxin(2:end)),0) || ~isequal(sum(naninxout(1:end-1)),0)
            error('Unexpected NaN in LickIn or LickOut.')  % NaN should be first LickIn or last LickOut
        end
        if wasnan && ~naninxin(1)
            error('Unexpected NaN in LickIn or LickOut.')  % NaN in LickOut should be followed by NaN in LickIn
        end
        if naninxout(end)
            wasnan = true;
        else
            wasnan = false;
        end
        TE.LickIn{k}(naninxin) = [];
        TE.LickOut{k}(naninxout) = [];
    else
        wasnan = false;
    end
end
if isempty(TE.LickIn{end})   % take care of different empty matrices which could result in errors later
    TE.LickIn{end} = zeros(0,1);
end
if isempty(TE.LickOut{end})
    TE.LickOut{end} = zeros(0,1);
end

% Save
if ifsave == 1
    
    % For each session, save the individual variables from the extracted TE
    savename = ['TE_' TE.Ratname{10} '_' TE.sessionID{2}];
    savename2 = 'TE';
    save([fileparts(filepath) filesep savename],'-struct','TE')
    save([fileparts(filepath) filesep savename2],'-struct','TE')
    
    % this save process saves individual variables in the file; can put into a
    % struct format by writing TE=load(TE_hr012_100718a.mat), and then can call
    % up TE.sessionID etc.
end
disp('Trial events structure completed.')