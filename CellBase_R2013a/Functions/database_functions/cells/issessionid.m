function yn = issessionid(sessionid)
%ISSESSIONID   Check if CellBase session ID exists.
%   YN = ISSESSIONID(SESSIONID) returns true if the referred session ID
%   exists and false if it doesn't. Session ID can be given in two format:
%   1-by-2 cell array containing animal and session IDs or a character
%   array containing the session ID without animal specification.
%
%   See also ISCELLID.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   27-Sept-2013

%   Edit log: BH 9/27/13, 5/20/14

% All sessions
allsessions = listtag('allsessions');

% Look for session ID
if iscell(sessionid)   % animalID and sessionID passed together
    inx1 = cellfun(@(s)isequal(s,sessionid{1,1}),allsessions(:,1));
    inx2 = cellfun(@(s)isequal(s,sessionid{1,2}),allsessions(:,2));
    yn = any(inx1&inx2);
else   % sessionID w/o animalID
    inx1 = cellfun(@(s)isequal(s,sessionid),allsessions(:,2));
    yn = any(inx1);
end