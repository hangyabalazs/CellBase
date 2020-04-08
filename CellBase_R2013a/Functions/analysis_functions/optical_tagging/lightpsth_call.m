function lightpsth_call(mousepath)
%LIGHTPSTH_CALL   Wrapper function for LIGHTPSTH.
%   LIGHTPSTH_CALL(MOUSEPATH) calls LIGHTPSTH for all sessions of an
%   animal.
%
%   See also LIGHTPSTH.

% Sessions
if nargin < 1
    mousepath = fullfile(getpref('cellbase','datapath'),'HDB18\');
end
dr = dir(mousepath);
dr = dr(3:end);

% Call 'lightpsth' for all sessions
NumSessions = length(dr);
H = nan(1,NumSessions);   % figure handles
for iS = 1:NumSessions
    lightpsth(fullfile(mousepath,dr(iS).name))
    H(iS) = gcf;
end
% keyboard