function [TRIALS, TAGS, RTAGS, NUM_TAGS, PartNum] = partition_trials(TE,partitions,varargin)
%PARTITION_TRIALS   Sort all trials to partitions.
%   [TRIALS, TAGS, RTAGS, NUM_TAGS, PARTNUM] = PARTITION_TRIALS(TE,PARTITIONS)
%   sorts all trials in TE structure to nonoverlapping groups according to
%   PARTITIONS. Partiotioned trials (TRIALS) are returned with
%   corresponding tags (TAGS, RTAGS) and number of elements in each
%   partition (PARTNUM).
%
%   [TRIALS, TAGS, RTAGS, NUM_TAGS, PARTNUM] = PARTITION_TRIALS(TE,PARTITIONS,VALIDTRIALS)
%   restricts the partitioning to valid trials.
%
%   See also VIEWCELL2B.

%   Edit log: BH 7/5/12, 9/4/12

% Input argument check
if nargin > 2
    validtrials = varargin{1};
else
    validtrials = [];
end

% Make partitions
allfields = fieldnames(TE);
if ischar(partitions)   % when there is only one
    partitions = {partitions};
end

iPART = 1;
lenpar = length(partitions);
PartNum = zeros(1,lenpar);
for iP = 1:lenpar
    PartNum(iP) = iPART - 1;   % number of partitions
    if strcmpi(partitions{iP},'all')
        TRIALS{iPART} = 1:length(TE.(allfields{1}));
        TAGS{iPART} = 'All';
        iPART = iPART + 1;
    else
        [exp2eval, values2match, condition, valid] = parse_token(partitions{iP},allfields);
        if strcmp(values2match,'none')
            if isempty(condition)
                TRIALS{iPART} = find(TE.(exp2eval));            
                TAGS{iPART} = exp2eval;
                RTAGS{iPART} = exp2eval;
            else
                TRIALS{iPART} = intersect(find(TE.(exp2eval)),find(TE.(condition)));  
                TAGS{iPART} = [exp2eval '&' condition];
                RTAGS{iPART} = [exp2eval '&' condition];
            end
            iPART = iPART + 1;
        else
            if isempty(values2match)
                avs = unique(TE.(exp2eval));
                avs = avs(~isnan(avs));
                values2match =  num2cell(avs);
                NUM_TAGS = [values2match{:}];
            else
                NUM_TAGS = [values2match{:}];
            end
            
            % Check if there are multiple options for conversion
            if iscell(values2match{1})
                allvals = unique(TE.(exp2eval));
                lenv2m = length(values2match);
                mismatches = zeros(1,lenv2m);
                for iL = 1:lenv2m
                    uvals = nonan(unique(cell2mat(values2match{iL})));
                    mismatches(iL) = length(setxor(allvals,uvals));
                end
                 values2tag = values2match{1};
                [junk, RULE2USE] = min(mismatches);
                values2match = values2match{RULE2USE};
            else
                values2tag = values2match;
            end
          
            for iL = 1:length(values2match)
                vals = mat2str(values2match{iL});
                newtag =  mat2str(values2tag{iL});  % always label it like the first description
                if isempty(condition)
                    TRIALS{iPART} = findpos(TE.(exp2eval),str2double(vals));
                    TAGS{iPART} = [exp2eval '=' newtag];
                    RTAGS{iPART} = [exp2eval '=' vals];
                else
                    TRIALS{iPART} = findpos(TE.(condition),1,TE.(exp2eval),str2double(vals));
                    TAGS{iPART} = [exp2eval '=' newtag '&' condition];
                    RTAGS{iPART} = [exp2eval '=' newtag '&' condition];
                end
                iPART = iPART + 1;
            end
        end
    end  
end % iP

% Restrict to valid trials
if ~isempty(validtrials)
    for iP = 1:length(partitions)
        TRIALS{iP} = intersect(TRIALS{iP},validtrials);
    end
end

% Number of elements in each partition ('PartNum'), tags
PartNum(iP+1) = iPART - 1; 
PartNum = diff(PartNum);
TAGS = noblank(TAGS);