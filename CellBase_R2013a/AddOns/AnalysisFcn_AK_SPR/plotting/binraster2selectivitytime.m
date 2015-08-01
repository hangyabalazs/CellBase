function D = binraster2selectivitytime(binraster,COMP,valid_trials,varargin);


% assign defaults if no arguments passed
default_args = { ...
        'Transform'             'swap';...
        'Cumulative'            'no'; ...
        'windowsize'            0.1; ...
        'shift'                 0.06; ...
        'dt'                    0.02; ...
        'FigureNum'             0; ...
    };

[g, error] = parse_args(default_args,varargin{:});


%positions{1}=selecttrial(TE,'Correct & WaterPokeValid')
%positions{2}=selecttrial(TE,'Error & WaterPokeValid')


%if valid_trials doesnt contain trial number then convert
if isbinary(valid_trials) 
     valid_trials = find(valid_trials);
end 


if iscell(COMP)
    for i = 1:length(COMP)
        positions{i} = intersect( COMP{i}, valid_trials);
    end %i
else
    error('huh');
end


if strcmpi(g.Cumulative,'yes')
    binraster=cumsum(binraster,2);
end

%-------------------------

W  = g.windowsize/g.dt;
dW = max(floor(g.shift/g.dt),1);

maxtpos = size(binraster,2);
LengthTime = ceil(maxtpos/dW);
for iD=1:LengthTime
    INC = (iD-1)*dW;
    br_sum(:,iD) = nansum(binraster(:,INC+1:min(W+INC,maxtpos))')';
end
br_sum(:,end) = nansum(binraster(:,end-5:end)')';
 
%---------------------

X=br_sum(positions{1},:);
Y=br_sum(positions{2},:);

Z=[X(:); Y(:)];

dmin = min(setdiff(diff(sort(Z)),0));
%nbin = ceil(max(length(x)*1.2,length(y)*1.2)); % some automatic assignment for number of bins
MN = min(Z);    MX = max(Z);
bin_size = dmin;
bins = MN-bin_size:bin_size:MX+bin_size;
Lx = size(X,1); Ly = size(Y,1);

MNtrials = 5;
for iT = 1:LengthTime  
    %D = auc(x,y,Lx,Ly,bins);
    %function D = auc(x,y,Lx,Ly,bins);
    x=X(:,iT);
    y=Y(:,iT);
    Lx = sum(~isnan(x)); Ly=sum(~isnan(y));
    if min([Lx Ly]) > MNtrials
        p = histc(x,bins);  q = histc(y,bins);
        cdf1 = cumsum(p)/Lx; cdf2 = cumsum(q)/Ly;
        D(iT)=trapz(cdf1,cdf2); 
    else
        D(iT) = NaN;
    end
end


%%%%%%%%%%%%%%
% transform D
%%%%%%%%%%%%%
switch lower(g.Transform)
    case {'swap','0.5','disc'}
        D=abs(D-0.5)+0.5;       % 'swap'
    case {'scale','selectivity'}
        D=2*(D-0.5);            % 'scale'
    otherwise
        D=D;            % no'ting to do;
end


if ~isempty(g.FigureNum) & g.FigureNum ~= 0
    figure(g.FigureNum)
    clf
    subplot(211)
    hold on;
    stairs(nanmean(binraster(positions{1},:)),'b');
    stairs(nanmean(binraster(positions{2},:)),'r'); 
    subplot(212)
    stairs(D,'k');
end

%
% return the right size
Ltime = size(binraster,2);
R=ceil(Ltime/length(D));
D2=interp(D,R);
D = D2(1:Ltime);
        
        