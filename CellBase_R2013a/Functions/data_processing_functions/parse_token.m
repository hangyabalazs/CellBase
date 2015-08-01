function [exp2eval, values2match, condition, valid] = parse_token(intoken,allfields)
%PARSE_TOKEN   Parse arguments.
%   [EXP2EVAL, VALUES2MATCH, CONDITION, VALID] = PARSE_TOKEN(INTOKEN,ALLFIELDS)
%   parses input arguments to output strings EXP2EVAL, VALUES2MATCH and
%   CONDITION. A validity check is performed with its result returned in
%   VALID.
%
%   Examples for INTOKEN:
%   intoken = '#OdorRatio:{[0 20 80 100] [32 44 56 69]} & Correct';
%   intoken = '#Stimulus & Correct';
%   intoken = '#OdorRatio:{{[0 100] [32 68] [44 56]} {[0 100] [0 100] [32 68]}}' 
%   intoken = '#OdorRatio:{[0 100] [32 44 56 68]}' 
%
%   See also PARSE_ARGS and PARTITION_TRIALS.

% Parse
intoken = deblank(intoken);
condition = '';
if intoken(1) == '#'
    [intoken2, condition] = strtok(intoken(2:end),'&');
    condition = strrep(strrep(condition,'&',''),' ','');
    
    [exp2eval, valuestring] = strtok(intoken2,':');
    if isempty(valuestring)
        values2match = [];
    else
        values2match = eval(valuestring(2:end));
    end
else
    exp2eval = intoken;
    values2match = 'none';
end

% Check validity
if exist('allfields','var')
    fieldname = strtok(exp2eval,' <>=~');
    valid = sum(ismember(allfields,fieldname));
else
    valid = NaN;
end