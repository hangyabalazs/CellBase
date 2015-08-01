function TE = solo2trialevents2_auditory_gonogo(filepath,ifsave)
%SOLO2TRIALEVENTS2_AUDITORY_GONOGO   Creates trial events structure.
%   SOLO2TRIALEVENTS2_AUDITORY_GONOGO(SESSPATH) creates trial events 
%   structure ('TrialEvents') for the session given in SESSPATH. The event
%   are extracted from the saved behavior file. The first trial is omitted.
%   Most of the events are given relative to the TrialStart event.
%
%   SOLO2TRIALEVENTS2_AUDITORY_GONOGO(SESSPATH,IFSAVE) uses an addidional
%   logical input argument do determine whether to save the results.
%
%   See also MAKETRIALEVENTS2_ANOGO and DEFINEEVENTSEPOCHS_AGONOGO.

%   Edit log: SPR 2011-03-03, BH 6/29/11, BH 7/7/11

% We want to save the result
if nargin < 2
    ifsave = 1;
end

% Determine filepath and load data
load(filepath)

% Session ID (date and 'a' or 'b') implement solofilename2tags here
sessionID = filepath(length(filepath)-10:length(filepath)-4);
sess_datetime = saved.SavingSection_SaveTime;

% No. of trials (exclude trial 1 and last trial may not have been completed
ntrials = saved.ProtocolsSection_n_completed_trials - 1;

% Preallocation of space
TE.TrialStart = nan(1,ntrials);
TE.TrialEnd = nan(1,ntrials);
TE.TrialsPerBlock = nan(1,ntrials);
TE.ChoiceSide = nan(1,ntrials);   % [0 1 2 --> no response, left, right]
TE.CorrectSide = nan(1,ntrials);   % [1 2 --> left, right (what about neither and both)]

TE.Hit = nan(1,ntrials);
TE.CorrectRejection = nan(1,ntrials);
TE.FalseAlarm = nan(1,ntrials);
TE.Miss = nan(1,ntrials);
TE.Omission = nan(1,ntrials);

TE.StimulusID = nan(1,ntrials);
TE.StimulusDuration = nan(1,ntrials);
TE.SoundIntensity = nan(1,ntrials);
TE.ResponseDelay = nan(1,ntrials);
TE.ResponseWindow = nan(1,ntrials);

TE.StimulusOn = nan(1,ntrials);
TE.StimulusOff = nan(1,ntrials);
TE.SignalModality = cell(1,ntrials);
TE.GoSignalModality = ones(1,ntrials) * 2;   % 0 no go signal, 1 LED, 2 Sound
TE.GoSignalOn = nan(1,ntrials);
TE.GoSignalOff = nan(1,ntrials);

TE.ReactionTime = nan(1,ntrials);
TE.TotalReactionTime = nan(1,ntrials);
TE.GoRT = nan(1,ntrials);
TE.NoGoRT = nan(1,ntrials);

TE.LeftPortIn = nan(1,ntrials);    % breaking the beam of the water port in headfixed protocol
TE.LeftWaterValveOn = nan(1,ntrials);
TE.LeftWaterVolume = nan(1,ntrials);
TE.LeftWaterValveDur = nan(1,ntrials);

TE.RightPortIn = nan(1,ntrials);
TE.RightWaterValveOn = nan(1,ntrials);
TE.RightWaterVolume = nan(1,ntrials);
TE.RightWaterValveDur = nan(1,ntrials);

TE.LightStimulation = nan(1,ntrials);
TE.LightStimulation2 = nan(1,ntrials);

TE.ErrorITI = nan(1,ntrials);
TE.ITIDistribution = nan(1,ntrials);
TE.ITIBegins = cell(1,ntrials);
TE.ITIEnds = cell(1,ntrials);
TE.TotalITI = nan(1,ntrials);
TE.NumITIRepeats = nan(1,ntrials);

TE.ResponseDelayStart = nan(1,ntrials);
TE.ResponseDelayEnd = nan(1,ntrials);
TE.ResponseWindowStart = nan(1,ntrials);
TE.ResponseWindowEnd = nan(1,ntrials);

TE.LeftPunishBegin = nan(1,ntrials);
TE.LeftPunishEnd = nan(1,ntrials);
TE.RightPunishBegin = nan(1,ntrials);
TE.RightPunishEnd = nan(1,ntrials);

TE.EarlyResponse = nan(1,ntrials);
TE.LeftEarlyResponse = nan(1,ntrials);
TE.RightEarlyResponse = nan(1,ntrials);
TE.Earlyresp = nan(1,ntrials);
TE.Earlyresp_inc = nan(1,ntrials);
TE.Earlyresp_2Light = nan(1,ntrials);
TE.Earlyresp_2Light_inc = nan(1,ntrials);
TE.Restart_itis = nan(1,ntrials);
TE.sessionID = cell(1,ntrials);
TE.datetime = nan(1,ntrials);

TE.LickIn = cell(1,ntrials);
TE.LickOut = cell(1,ntrials);

% Correct side
lefttrials = find(saved.auditory_gonogo_TrialTypeList(1:ntrials)=='l');
righttrials = find(saved.auditory_gonogo_TrialTypeList(1:ntrials)=='r');
TE.CorrectSide(lefttrials) = 1; %#ok<FNDSB>
TE.CorrectSide(righttrials) = 2; %#ok<FNDSB>
TE.CorrectSide = [TE.CorrectSide(2:end) NaN];  % remove first trial

% Type of trial (for e.g. diff stim durations get different index)
TE.StimulusID = saved.auditory_gonogo_TrialTypes(1:ntrials);
TE.StimulusID = [TE.StimulusID(2:end) NaN];  % remove first trial

% Stimulus duration / sound intensity on every trial
TE.StimulusDuration = saved.auditory_gonogo_StimDurList(1:ntrials);
TE.StimulusDuration = [TE.StimulusDuration(2:end) NaN];  % remove first trial
TE.SoundIntensity = saved.auditory_gonogo_StimDurList(1:ntrials);   % it is actually sound intensity; however, I keep the redundant StimulusDuration field to avoid errors in other programs - BH
TE.SoundIntensity = [TE.SoundIntensity(2:end) NaN];  % remove first trial

% Animal name
TE.Ratname = saved_history.SavingSection_ratname(1:ntrials)';
TE.Ratname = [TE.Ratname(2:end) NaN];  % remove first trial

try
    % Light stimulation on each trial
    TE.LightStimulation = saved.auditory_gonogo_LightStimList(1:ntrials);
    TE.LightStimulation = [TE.LightStimulation(2:end) NaN];  % remove first trial
    TE.LightStimulation2 = saved.auditory_gonogo_LightStimulation2(1:ntrials);
    TE.LightStimulation2 = [TE.LightStimulation2(2:end) NaN];  % remove first trial
catch
end

% Define variables trialwise
countr = 0;
for iT = 2:ntrials
    countr = countr + 1;
    try TE.StimulusOn(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.deliverstim(1);catch end
    try TE.StimulusOff(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.deliverstim(2);catch end
    try TE.GoSignalOn(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.deliverstim(1);catch end
    try TE.GoSignalOff(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.deliverstim(2);catch end
%     try TE.ReactionTime(countr) = diff(saved_history.ProtocolsSection_parsed_events{iT}.states.beginresponse);catch end
    
    % Hit and Correct Rejection
    if ~isempty(saved_history.ProtocolsSection_parsed_events{iT}.states.correctreject)
        TE.CorrectRejection(countr) = 1;
        TE.ChoiceSide(countr) = 0;  % correct nogo response made
    end
    if ~isempty(saved_history.ProtocolsSection_parsed_events{iT}.states.hit)
        TE.Hit(countr) = 1;   % correct Go response
        TE.ChoiceSide(countr)=1;    % left choice made
        try
            TE.LeftPortIn(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.hit(1);
            TE.TotalReactionTime(countr) = TE.LeftPortIn(countr) - TE.StimulusOn(countr);
            TE.ReactionTime(countr) = TE.LeftPortIn(countr) - TE.GoSignalOn(countr);
            TE.GoRT(countr) = TE.ReactionTime(countr);
        catch
        end    % identical to pokes.L(1)
    end
    
    % False Alarm
    if ~isempty(saved_history.ProtocolsSection_parsed_events{iT}.states.falsealarm)
        TE.FalseAlarm(countr) = 1;
        TE.ChoiceSide(countr) = 1;
        try 
            TE.LeftPortIn(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.falsealarm(1);
            TE.TotalReactionTime(countr) = TE.LeftPortIn(countr) - TE.StimulusOn(countr);
            TE.ReactionTime(countr) = TE.LeftPortIn(countr) - TE.GoSignalOn(countr);
            TE.NoGoRT(countr) = TE.ReactionTime(countr);            
        catch
        end    % identical to pokes.L(1)
    end
    
    if ~isempty(saved_history.ProtocolsSection_parsed_events{iT}.states.miss)
        TE.Miss(countr) = 1;
        TE.ChoiceSide(countr) = 0;
    end
    
    % Trial Start (state 0)
    try TE.TrialStart(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.state_0(1,2); catch end
    % Trial End (ending state 0)
    try TE.TrialEnd(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.state_0(2,1); catch end
    try TE.ITIDistribution(countr) = saved.auditory_gonogo_ITIs(iT); catch end
    try TE.ITIBegins{countr} = saved_history.ProtocolsSection_parsed_events{iT}.states.iti(:,1); catch end
    try TE.ITIEnds{countr} = saved_history.ProtocolsSection_parsed_events{iT}.states.iti(:,2); catch end
    try TE.TotalITI(countr) = TE.ITIEnds{countr}(end) - TE.ITIBegins{countr}(1); catch end
    try TE.NumITIRepeats(countr) = length(TE.ITIBegins{countr}); catch end
    
    % Stimulus
    try TE.SignalModality{countr} = saved.auditory_gonogo_SignalModality; catch end
    try TE.StimulusOn(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.deliver_stim1(1); catch end
    try TE.StimulusOff(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.deliver_stim1(2); catch end
    
    try TE.ResponseDelayStart(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.waitforresponse(1); catch end
    try TE.ResponseDelayEnd(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.waitforresponse(2); catch end
    try TE.ResponseWindowStart(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.beginresponse(1); catch end
    try TE.ResponseWindowEnd(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.beginresponse(2); catch end
    try TE.GoSignalOn(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.beginrespsound(1); catch end
    try TE.GoSignalOff(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.beginrespsound(2); catch end
    
    try TE.LeftWaterValveOn(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.ldeliver(1); catch end
    try TE.LeftWaterVolume(countr) = saved_history.WaterValvesSection_Left_volume{iT}; catch end
    try TE.LeftWaterValveDur(countr) = saved_history.WaterValvesSection_LeftWValveTime{iT}; catch end

    try TE.RightWaterValveOn(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.rdeliver(1); catch end
    try TE.RightWaterVolume(countr) = saved_history.WaterValvesSection_Right_volume{iT}; catch end
    try TE.RightWaterValveDur(countr) = saved_history.WaterValvesSection_RightWValveTime{iT}; catch end

    try TE.LeftPunishBegin(countr) = saved_history.ProtocolsSection_parsed_events{2}.states.lpunish(1); catch end
    try TE.LeftPunishEnd(countr) = saved_history.ProtocolsSection_parsed_events{2}.states.lpunish(2); catch end
    try TE.RightPunishBegin(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.rpunish(1); catch end
    try TE.RightPunishEnd(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.rpunish(2); catch end
    
    try TE.EarlyResponse(countr) = saved_history.ProtocolsSection_parsed_events{iT}.states.earlyresponse(1); catch end
    try TE.sessionID{countr} = sessionID; catch end
    try TE.LickIn{countr} = saved_history.ProtocolsSection_parsed_events{iT}.pokes.L(:,1); catch end
    try TE.LickOut{countr} = saved_history.ProtocolsSection_parsed_events{iT}.pokes.L(:,2); catch end
end

% Sound duration
TE.SoundDuration = nan(1,ntrials);
TE.SoundDuration(:) = saved.auditory_gonogo_SoundDuration;
TE.SoundDuration = [TE.SoundDuration(2:end) NaN];  % remove first trial

% Sound Frequency
TE.SoundFrequency = nan(1,ntrials);
TE.SoundFrequency(:) = saved.auditory_gonogo_SoundFrequency;
TE.SoundFrequency = [TE.SoundFrequency(2:end) NaN];  % remove first trial

% Response
TE.ResponseDelay = nan(1,ntrials);
TE.ResponseDelay(:) = saved.auditory_gonogo_Delay2Resp;
TE.ResponseDelay = [TE.ResponseDelay(2:end) NaN];  % remove first trial

TE.ResponseType = nan(1,ntrials);
TE.ResponseType(TE.Hit==1) = 1;    % because TE.Hit is already OK, we don't have to remove first trial here
TE.ResponseType(TE.FalseAlarm==1) = 2;
TE.ResponseType(TE.CorrectRejection==1) = 3;
TE.ResponseType(TE.Miss==1) = 4;

% Assign a Block Number
% we are going to assume there was only 1 Block size ... .i.e. only one value for
% TrialsPerBlock per session.
TE.TrialsPerBlock(:) = saved_history.auditory_gonogo_TrialsPerBlock{1};
TE.TrialsPerBlock = [TE.TrialsPerBlock(2:end) NaN];  % remove first trial

blocklength = saved_history.auditory_gonogo_TrialsPerBlock{1}*length(unique(TE.StimulusDuration))*2;
numblocks = ceil(ntrials/blocklength);   % each block has equal number of A+ and B- trials of each type
blockid = reshape(repmat(1:numblocks,blocklength,1),1,[]);
TE.BlockNum = blockid(1:ntrials);
TE.BlockNum = [TE.BlockNum(2:end) NaN];  % remove first trial

% Convert to relative time to TrialStart
TE.StimulusOn = TE.StimulusOn - TE.TrialStart;
TE.StimulusOff = TE.StimulusOff - TE.TrialStart;
TE.ResponseDelayStart = TE.ResponseDelayStart - TE.TrialStart;
TE.ResponseDelayEnd = TE.ResponseDelayEnd - TE.TrialStart;
TE.ResponseWindowStart = TE.ResponseWindowStart - TE.TrialStart;
TE.ResponseWindowEnd = TE.ResponseWindowEnd - TE.TrialStart;
TE.GoSignalOn = TE.GoSignalOn - TE.TrialStart;
TE.GoSignalOff = TE.GoSignalOff - TE.TrialStart;
TE.LeftPortIn = TE.LeftPortIn - TE.TrialStart;
TE.LeftWaterValveOn = TE.LeftWaterValveOn - TE.TrialStart;
TE.RightPortIn = TE.RightPortIn - TE.TrialStart;
TE.RightWaterValveOn = TE.RightWaterValveOn - TE.TrialStart;

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
disp('TE analysis complete')