function out = meanrate(cellid,varargin)
%
%  MEANRATE     Calculates the mean firing of a cell
%
%
%   x = meanrate(cellid,{isi_threshold})
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Analysis functions should always take the {cellid} as
%  the first argument and can possibly take a variable argument
%  list. 
%  The output can be a *single* scalar or a vector of numerical values.
%  (i.e. don't give multiple outputs, put it into a vector)
%
%  Calling stimes = loadcell(cellid,'caller'); initially will load 
%  the necessary data.
%  
%  To support self-describing analysis functions return a cell array
%  of property descriptor strings back if cellid is 'default'.
%

   if strcmpi(cellid,'default')
       out = 'Rate';
       % If you define default property names make sure there are as many
       % as regular outputs. Use cellstr arrays for multiple outputs:
       % out = {'Rate','Mean ISI'}
       return
   end
   %
   stimes = loadcb(cellid,'Spikes');
    %
    %%
    %%%
    isi=diff(stimes);
    if nargin == 1
        out = 1/mean(isi);   % rate
    elseif nargin == 2 
        % the argument is always in the first varargin & could be a vector
        % or a cell array, so make sure you take out only the relevant
        % things
        argument = double(varargin{1});
        isi_thres = argument(1); %in milliseconds now
        out = 1/mean(isi(find(isi>isi_thres)));    
    end
