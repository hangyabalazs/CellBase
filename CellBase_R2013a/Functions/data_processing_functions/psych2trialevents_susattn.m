function TE = psych2trialevents_susattn(filepath,ifsave)
%PSYCH2TRIALEVENTS_SUSATTN   Creates trial events structure.
%   PSYCH2TRIALEVENTS_SUSATTN(SESSIONPATH) creates TrialEvents structure
%   for the human sustained attention session passed in SESSIONPATH. The
%   events are extracted from the saved behavior file. All time stamps are
%   stored relative to the TrialStart event.
%
%   PSYCH2TRIALEVENTS_SUSATTN(SESSIONPATH,IFSAVE) uses an addidional
%   logical input argument do determine whether to save the results.
%
%   See also SUSATTN_MAIN.

%   Edit log: BH 7/31/12

% We want to save the result
if nargin < 2
    ifsave = 1;
end

% Determine filepath and load data
load(filepath)

% Session ID
sessionID = filepath(length(filepath)-10:length(filepath)-4);
[pathname filename] = fileparts(filepath); %#ok<ASGLU>
cmps = strread(filename,'%s','delimiter','_');
subject = cmps{1};

% No. of trials (exclude trial 1 and last trial may not have been completed
[nblocks ntrialsperblock] = size(history.trial.TrialStart);
ntrials = nblocks * ntrialsperblock;
ntrials_completed = sum(~isnan(history.trial.TrialStart(:)));

% General
TE = struct;   % preallocate
TE.SubjectName = repmat({subject},1,ntrials);   % subject initials
TE.Ratname = TE.SubjectName;    % maintain compatibility with CellBase
TE.sessionID = repmat({sessionID},1,ntrials);   % session ID
if isfield(history.session,'SessionStart')    % session start time
    TE.datetime = repmat({history.session.SessionStart},1,ntrials);
end
TE.BlockNum = reshape(repmat(1:nblocks,ntrialsperblock,1),1,[]);

% Trial timing
TE.TrialStart = reshape(history.trial.TrialStart',1,ntrials);   % trial start
TE.TrialEnd = [TE.TrialStart(2:end) NaN] - TE.TrialStart;   % trial end, relative to trial start

% Fix trial variables
TE.TrialsPerBlock = repmat(ntrialsperblock,1,ntrials);  % number of trials in each block
if isfield(history.session,'NoPressITI')
    TE.NoPressITI = repmat(history.session.NoPressITI,1,ntrials);   % length of NoPressITI state (box on)
else
    TE.NoPressITI = repmat(1.5,1,ntrials);   % not recorded for the first session
    warning('psych2trialevents_susattn:NOPRESSITI','No record for NoPressITI.')
end

% Outcome variables
TE.Hit = zero2nan(reshape(history.trial.Hit',1,ntrials));   % HIT
TE.CorrectRejection = zero2nan(reshape(history.trial.CorrectReject',1,ntrials));  % CORRECT REJECTION
TE.FalseAlarm = zero2nan(reshape(history.trial.FalseAlarm',1,ntrials));   % FALSE ALARM
TE.Miss = zero2nan(reshape(history.trial.Miss',1,ntrials));   % MISS
TE.Early = zero2nan(reshape(history.trial.Early',1,ntrials));   % restart ITIs (new trial)
TE.ResponseType = nan(1,ntrials);   % use integer codes for the four possible outcomes
TE.ResponseType(TE.Hit==1) = 1;   % Hit = 1
TE.ResponseType(TE.FalseAlarm==1) = 2;   % FA = 2
TE.ResponseType(TE.CorrectRejection==1) = 3;  % CR = 3
TE.ResponseType(TE.Miss==1) = 4;   % Miss = 4
TE.ResponseType(TE.Early==1) = 5;   % Restart ITI = 5

% Reaction time
TE.ReactionTime = reshape(history.trial.RT',1,ntrials);   % reaction time as recorded in the task
TE.GoRT = TE.ReactionTime;   % RT for Hit trials only (NaN for other trials)
TE.GoRT(isnan(TE.Hit)) = NaN;
TE.NoGoRT = TE.ReactionTime;
TE.NoGoRT(isnan(TE.FalseAlarm)) = NaN;   % 'RT' for FA trials only (NaN for other trials)
TE.LeftPortIn = reshape(history.trial.ResponseTime',1,ntrials);  % 'response time' relative to trial start

% Foreperiod variables
if size(history.trial.ITI,1) > 1
    ITI = history.trial.ITI;  % foreperiod length
else
    ITI = repmat(history.trial.ITI,nblocks,1);   % in the initial version, ITIs were repeated for the blocks
end
TE.NoPokeITIBegins = zeros(1,ntrials);   % start of NoPressITI period, box on (same as trial start, time 0 rel. to trial start)
TE.NoPokeITIEnds = TE.NoPokeITIBegins + TE.NoPressITI;   % end of NoPressITI period, box off
TE.ITIDistribution = reshape(ITI',1,ntrials);    % foreperiod length
TE.ITIBegins = TE.NoPokeITIEnds;   % start of foreperiod, box off
TE.ITIEnds = TE.ITIBegins + TE.ITIDistribution;   % end of foreperiod (should equal stimulus on)

% Stimulus variables
TE.StimulusID = reshape(history.trial.TrialType',1,ntrials);  % 1 = go trial; 0 = no-go trial
TE.SoundDuration = repmat(history.session.StimDur,1,ntrials);  % preset tone duration (but note: terminated by response - see StimulusDuration)
TE.SoundFrequency = nan(1,ntrials);   % frequency of the stimulus tone
TE.SoundFrequency(TE.StimulusID==1) = history.session.GoFreq;   % go tone frequency
TE.SoundFrequency(TE.StimulusID==0) = history.session.NoGoFreq;   % no-go tone frequency
TE.StimulusDuration = min(repmat(history.session.StimDur,1,ntrials),...
    TE.ReactionTime);   % tone duration; response terminates the tone
TE.StimulusDuration(TE.Early==1) = NaN;   % no stimulus was played in Early trials
TE.SoundIntensity = reshape(history.trial.StimAmplitude',1,ntrials);   % tone intensity
if isfield(history.trial,'SimulusOn')
    TE.StimulusOn = reshape(history.trial.StimulusOn',1,ntrials);  % tone onset time relative to trial start
else
    TE.StimulusOn = TE.ITIEnds;   % initially no time stamp was recorded for tone onset
end
TE.StimulusOff = TE.StimulusOn + TE.StimulusDuration;   % tone offset
TE.SignalModality = repmat({'AUDITORY'},1,ntrials);   % signal modality = 'AUDITORY'
TE.GoSignalModality = repmat(2,1,ntrials);  % 2 = Sound (1 = LED)

% Response window
TE.ResponseWindow = repmat(history.session.StimDur,1,ntrials) +...
    repmat(history.session.ResponseWait,1,ntrials);   % a response is allowed durung the tone and a short period after it, termed 'ResponseWait'
TE.ResponseWindowStart = TE.StimulusOn;   % start of response window: tone onset
TE.ResponseWindowEnd = TE.ResponseWindowStart + TE.ResponseWindow;  % end of response window

% Include only completed trials
fld = fieldnames(TE);  % field names of the trial events structure
nfld = length(fld);   % number of fields
for k = 1:nfld   % loop through the fields
    TE.(fld{k}) = TE.(fld{k})(1:ntrials_completed);  % restrict to completed trials
end

% Save
if ifsave == 1
    
    % Save the individual variables for each session 
    savename = ['TE_' TE.SubjectName{2} '_' TE.sessionID{2}];
    savename2 = 'TE';   % save under two different names to maintain compatibility with all CellBase versions
    save([fileparts(filepath) filesep savename],'-struct','TE')
    save([fileparts(filepath) filesep savename2],'-struct','TE')
    
    % this saves fields of the trial events structure as individual
    % variables in the file; reconstruct structure format by calling TE =
    % load(FILENAME) (e.g. TE = load(TE_bh_120727h.mat) for loading trial
    % events
end
disp('TE analysis complete')