function  print_progress(t,TPRINT,N)
%PRINT_PROGRESS   Progress indicator.
%   PRINTPROGRESS(T,TPRINT,N) indicates progress in Command Window. T is
%   the serial number of computation step, TPRINT controls the frequency of
%   Command Window displays and N controls the frequency of numeric
%   displays (other displays are single dots).
%
%   Example:
%   print_progress(i,round(n/100),5)
%
%   See also WAITBAR.

%   Edit log: BH 3/23/11

% Print to Command Window
if ~mod(t,TPRINT)
    if ~mod(t,TPRINT*N)
        fprintf('%d',t/TPRINT)
    else
        fprintf('.')
    end
end