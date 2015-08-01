function [speedpsth, speedpsth_sd, speedspsth, speedspsth_se] = makespeedpsth(speed,postime,evtimes,time,margins,dt,sigma,COMP,valid_trials)
%[SPEEDPSTH, SPEEDSPSTH, SPEEDSPSTH_SE] = SPEEDPSTH(SPEED,EVTIMES,TIME,MARGIN,G.DT,G.SIGMA,COMPTRIALS,VALID_TRIALS)
% makes a psth of mouse speed aligned to an event
% Can substitute speed with angle, X or Y position.
% SPR 2011-01-05

% Get event aligned speed psth matrix for all trials.
NumTrials = length(evtimes);
speedpsthmatrix = nan(NumTrials,length(time));
for iT = 1:length(evtimes),
    tlimit = [evtimes(iT)+time(1) evtimes(iT)+time(end)];
    ind = find(postime >= tlimit(1) & postime <= tlimit(2));
    X = postime(ind);
    Y = speed(ind);
    Xi = tlimit(1):dt:tlimit(2);
    Yi = interp1(X,Y,Xi,'linear','extrap');
    speedpsthmatrix(iT,:) = Yi;
end

% Get trial indices for each partition.
if iscell(COMP)
    nCOND = length(COMP);
    for i = 1:nCOND
        positions{i} = intersect( COMP{i}, valid_trials);
    end %i
else % matrix
    iCOND = 1;
    nCOND = size(COMP,1);
    for i = 1:nCOND
        CONDITIONS = unique(COMP(i,:));
        if length(CONDITIONS) > 2 && ~ismember(100,CONDITIONS)
            CONDITIONS = setdiff(CONDITIONS, 0);  %conditions to compare, 0 doesn't count
        end
        for j = 1:length(CONDITIONS)
            positions{iCOND} = intersect( find(COMP(i,:) == CONDITIONS(j)), valid_trials);
            iCOND = iCOND + 1;
        end %j
    end %i
end %if

speedpsth     = nan(nCOND,size(speedpsthmatrix,2));
speedpsth_sd  = speedpsth;
speedspsth    = speedpsth;
speedspsth_se = speedpsth;

for iCOND = 1:length(positions)
    
    NUM_TRIALS = length(positions{iCOND});  % doesn't account for variable windows
    
    if ~isempty(positions{iCOND})  
        speedpsth(iCOND,:)     = nanmean(speedpsthmatrix(positions{iCOND},:)); 
        speedpsth_sd(iCOND,:)  = nanstd(speedpsthmatrix(positions{iCOND},:));
        speedspsth(iCOND,:)    = smoothed_psth(speedpsthmatrix(iCOND,:),dt,sigma);
        speedspsth_se(iCOND,:) = smoothed_psth(speedpsth_sd(iCOND,:)./sqrt(NUM_TRIALS-1),dt,sigma);
    end
end %iCOND