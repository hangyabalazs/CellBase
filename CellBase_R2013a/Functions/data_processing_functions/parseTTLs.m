function [ttls2 onttl offttl] = parseTTLs(ttls)
%PARSETTLS   Parse TTLs to channels.
%   [TTLS2 ONTTL OFFTTL] = PARSETTLS(TTLS) converts TTLs to binary numbers
%   (TTLS2) and returns TTLs onsets (ONTTL) and offsets (OFFTTL) of
%   different channels.
%
%   See also MAKETRIALEVENTS2_GONOGO.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   25-Sep-2012

%   Edit log: BH, 9/25/12

% Convert to binary
ttls2 = dec2bin(ttls);

% Future and past TTL
[n dpth] = size(ttls2);   % number of TTLs and bit depth
past_ttl = [zeros(1,dpth); ttls2(1:end-1,:)];  % series of past TTLs
future_ttl = [ttls2(2:end,:); zeros(1,dpth)];  % series of future TTLs

% TTL onset and offset
onttl = nan(n,dpth);
offttl = nan(n,dpth);
for k = 1:size(ttls2,2)
    onttl(:,k) = past_ttl(:,k)==0 & ttls2(:,k)==1;   % TTL onset
    offttl(:,k) = past_ttl(:,k)==1 & ttls2(:,k)==0;   % TTL offset
end

% -------------------------------------------------------------------------
function s = dec2bin(d)

d = d(:);
[f e] = log2(max(d));
s = rem(floor(d*pow2(1-e:0)),2);   % convert to binary