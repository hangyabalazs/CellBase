function poscommon = findpos(varargin)
%FINPOS   Find trials with specific properties.
%   POSCOMMON = FINDPOS(TRIALS,VALUE) finds the indices of trials
%   (POSCOMMON) which have matching properties specified by the TRIALS,
%   VALUE pairs. Any number of such pairs can be passed. TRIALS should
%   contain the field of the trial events structure which is to be
%   searched for the specific values. 
%
%   Example:
%   pos = findpos(TE.ResponseType,[1 2],TE.StimDuration,40)
%
%   See also PARTITION_TRIALS.

%   Edit log: BH 6/23/11

% Input argument check
if mod(nargin,2)
    disp('Incorrect args')
end

% Find positions
for iARG = 0:(nargin/2-1);
    X = varargin{iARG*2+1};
    TOFIND = varargin{iARG*2+2};
    posfound{iARG+1} = [];
    for iF = 1:length(TOFIND)
        posnew = find(X == TOFIND(iF));
        posfound{iARG+1} = union(posfound{iARG+1},posnew);
    end
end

% Find intersect
poscommon = 1:length(X);
for iCOND = 1:length(posfound);
    poscommon = intersect(poscommon,posfound{iCOND});
end