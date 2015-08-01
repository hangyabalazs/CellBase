function sessionids = findsession(varargin)
%FINDSESSION   Locate sessions in CellBase.
%   FINDSESSION returns the session ID(s) in CellBase for a particular
%   animal or date. A list of animals (cell array of animal IDs) is also
%   supported. Date should be passed in YYMMDD format. The session IDs are
%   returned along with the corresponding animal IDs.
%
%   Syntax:
%   SESSIONIDS = FINDSESSION('ANIMAL',ANIMALID)
%   SESSIONIDS = FINDSESSION('DATE',DATE)
%   SESSIONIDS = FINDSESSION('ANIMAL',ANIMALID,'DATE',DATE)
%
%   See also FINDCELL.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   27-Sept-2013

%   Edit log: BH 9/27/13

% Input arguments
prs = inputParser;
addParamValue(prs,'animal',{''},@(s)ischar(s)|iscellstr(s))   % animal ID
addParamValue(prs,'rat',{''},@(s)ischar(s)|iscellstr(s))   % rat ID
addParamValue(prs,'mouse',{''},@(s)ischar(s)|iscellstr(s))   % mouse ID
addParamValue(prs,'date',{''},@(s)ischar(s)|iscellstr(s))   % date
parse(prs,varargin{:})
g = prs.Results;
fld = fieldnames(g);
cellfun(@convertfield,fld);  % convert all fields to cells
    function convertfield(flnm)   % convert to cell
        if ~iscell(g.(flnm))
            g.(flnm) = {g.(flnm)};
        end
    end

% Load CellBase
global CELLIDLIST ANALYSES TheMatrix
if isempty(CELLIDLIST)
    load(getpref('cellbase','fname'));
end

% All sessions
sessionids = listtag('sessions');

% Select animals
if ~isempty(g.animal{1}) || ~isempty(g.mouse{1}) || ~isempty(g.rat{1})
    inx = cellfun(@(s)ismember(s,[g.animal g.mouse g.rat]),sessionids(:,1));
    sessionids = sessionids(inx,:);
end

% Select date
if ~isempty(g.date{1})
    inx = cellfun(@(s)ismember(s(1:end-1),g.date),sessionids(:,2));
    sessionids = sessionids(inx,:);
end
end