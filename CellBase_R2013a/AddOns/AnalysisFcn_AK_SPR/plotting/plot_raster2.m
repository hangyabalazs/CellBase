function fhandle0 = plot_raster2(stimes,time,valid_trials,COMPTRIALS,TAGS,EventTimes,window_margin,ev_windows,sort_var,params,varargin)
% This can plot a raster on a specified axis/subplot within a figure

%plot_raster(fighandle,time,binraster,trial_order,EventTimes,tlimits,partitions,partition_colors,labelx,labely);


%fill([trig trig epoch epoch],ydata,shade_color,'Edgecolor',shade_color);

%spikes = ST.event_stimes{9};

%stimes2binraster(stimes,time,dt,ev_windows,window_margin);


default_args = params;
default_args.TitleStr = '';
default_args.XLabel = 'Time';
default_args.PSTH = 'on';
default_args.NumTrials2Plot = 15;
default_args.RasterIDbarDX  = 0.01;
default_args.RasterIDbarWidth = 6;
default_args.RasterPanelDistance = 0.08;
default_args.RasterPlotNumtrials = 1;
default_args.EventMarkerWidth=0.008;
default_args.EventMarkerHeight=2;
default_args.SpikeWidth = 0.5;
default_args.ForceSort=1; % sort even in the presence of NaN values (good for Previous and NExt Events where either the first or last value is a NaN

[par, error] = parse_args(default_args,varargin{:});



NumTrials = length(stimes);
if ~exist('window_margin','var')
   window_margin = [0 0];
end
if exist('ev_windows','var')
    if size(ev_windows,1) ~= NumTrials
        ev_windows = ev_windows';
    end
    ev_windows = ev_windows + repmat(window_margin,NumTrials,1);
end



det = par.EventMarkerWidth/2;
dnt = par.EventMarkerHeight/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumParts = length(COMPTRIALS);

if strcmpi(par.PSTH,'on')
    NumPanels = NumParts + 1;
else
    NumPanels = NumParts;
end

% if strcmp(get(par.FigureNum,'Type'),'axes'),
if ishandle(par.FigureNum),
    if strcmp(get(par.FigureNum,'Type'),'figure')
        clf;
        fhandle0=axes;
    elseif strcmp(get(par.FigureNum,'Type'),'axes')
        fhandle0=par.FigureNum;
        cla;
    end
    ax=splitaxes(fhandle0,NumPanels);
else
    fhandle0=focusfigure(par.FigureNum,'create');
    fhandle0=axes;
    ax=splitaxes(fhandle0,NumPanels);
%     ax=subplot(NumPanels,1,NumPanels);
end

iAX=1;

for iCOMP=1:NumParts
    
    %get trials to plot
    ind = intersect(COMPTRIALS{iCOMP},valid_trials);    %use only valid ones
    rndind = randperm(length(ind));
    ind = ind(rndind(1:min(length(ind),par.NumTrials2Plot)));
    NumTrials(iCOMP) = length(ind);
        %reindex
        if par.ForceSort==1,
            [i,sind]=sort(sort_var(ind));                       %sort by the sort variable
            ind = ind(sind);
        else % do it only if there are no NaN values
            if ~isnan(sort_var) %reindex
                [i,sind]=sort(sort_var(ind));                       %sort by the sort variable
                ind = ind(sind);
            end
        end %  
    
    if NumTrials(iCOMP) > 0
        
        nt=1:NumTrials(iCOMP);
%         subplot(NumPanels,1,iCOMP)
        subplot(ax(iAX));
        hold on
        
        for iTRIAL=1:NumTrials(iCOMP)
            indTRIAL=ind(iTRIAL);
            spikes = stimes{indTRIAL}; % SPR 2009/12/24: The stimes format is in the right way. No need to transpose
            win = [max(ev_windows(indTRIAL,1),time(1)) min(ev_windows(indTRIAL,2),time(end))]; 
            ok_spikes = spikes(find(spikes > win(1) & spikes <= win(2)));
            %ind_ok_spikes = ceil((ok_spikes-time(1))/dt);       %why subtract? 
            if ~isempty(ok_spikes)
                if size(ok_spikes,1)==1,
                    line([ok_spikes; ok_spikes],[(nt(iTRIAL)-0.5)*ones(1,length(ok_spikes)); (nt(iTRIAL)+0.5)*ones(1,length(ok_spikes))],'Color','k','Linewidth',par.SpikeWidth);
                else
                    line([ok_spikes ok_spikes]',[(nt(iTRIAL)-0.5)*ones(1,length(ok_spikes)); (nt(iTRIAL)+0.5)*ones(1,length(ok_spikes))],'Color','k','Linewidth',par.SpikeWidth);
                end
            end
        end %iTrial
        
        for iE = 1:size(EventTimes,1)
            %plot(EventTimes(iE,:),1:NumTrials, gcolor(iE+5,'.'))
            et=EventTimes(iE,ind);
            et(find(et<time(1)))   = NaN;
            et(find(et>time(end))) = NaN;
            %            
            %if (length(find(et==0)) ==  length(nonan(et)))  %All the not-NaN values are zeros
            if length(find(et==0)) == NumTrials(iCOMP)
                line([0 0],[1 NumTrials(iCOMP)],'Color',par.ShowEventsColors{iE},'LineWidth',1);
            else
                p(iE)=patch([et-det; et-det; et+det; et+det],[nt+dnt; nt-dnt; nt-dnt; nt+dnt],'k');
                try
                    set(p(iE),'FaceColor',par.ShowEventsColors{iE},'EdgeColor','none');
                catch
                    set(p(iE),'FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',0.1);
                end
            end
        end %iE
%         axis([time(1) time(end) 1-2*eps NumTrials(iCOMP)]);
        axis([par.window(1) par.window(2) 1-2*eps NumTrials(iCOMP)]);
%         ax(iAX)=gca; 
        if (par.RasterPlotNumtrials)
            set(ax(iAX),'YAxisLocation','left','YColor','k','YTick',[1-2*eps NumTrials(iCOMP)],'YTickLabel',{'',num2str(NumTrials(iCOMP))});
        end
        if NumTrials(iCOMP) > 1
            set(ax(iAX),'XTickLabel','','XColor','w');
            set(ax(iAX),'YAxisLocation','left','XColor','w');           
        end
         yl(iAX) = ylabel(TAGS{iCOMP});
        iAX=iAX+1;
    end % NumTrials > 0
end %iCOMP

%-------------------------------------------------
NewIndex = find(NumTrials~=0);
NumTrials = NumTrials(NewIndex);
NumParts = length(NumTrials);

%--------------------------------------------------------------------------
%Put up panel for PSTH
if strcmpi(par.PSTH,'on')
    fhandle0(2)=subplot(NumPanels,1,NumPanels);
    ax(iAX) = gca;
    NumPanels = NumParts+1;
else
    set(ax(end),'XTickLabelMode','auto','XColor','k','FontName','Ariel','FontSize',12);
end %iCOMP

%set(ax,'TickDir','out','YTick',[1 par.NumTrials2Plot],'YTickLabel','');
set(ax,'TickDir','out');
set(yl,'FontName','Ariel','FontSize',10);

%-------------------------------------------------
%REARRANGE
%
if ~isnan(par.NumTrials2Plot)
    
    %Rearrange equally
    pos=get(ax(1),'Position');
    for i=2:NumPanels
        set(ax(i),'Position',pos-[0 (i-1)*(pos(4)*(1+par.RasterPanelDistance)) 0 0]);
    end
    
else   
    posS=get(ax(1),'Pos');
    posE=get(ax(NumParts),'Pos');
    
    %Distance between panels
    deltaPD=ceil(posS(4)*par.RasterPanelDistance*1000)/1000;
    %deltaY = posS(2)+posS(4)-posE(2)-(NumParts-1)*deltaPD-0.1;
    deltaY = 0.6;
    
    TotalTrials=sum(NumTrials);
    %pbottom(NumParts)=posE(2)-0.05;
    pbottom(NumParts) = 0.06+0.2+deltaPD;
    for i=NumParts:-1:1
        pheight(i) = max(ceil(deltaY*NumTrials(i)/TotalTrials*100)/100,deltaPD);
        set(ax(i),'Position',[posS(1) pbottom(i) posS(3) pheight(i)]);
        if i>1
            pbottom(i-1) = pbottom(i) + pheight(i) + deltaPD;
        end
    end %i
end %NumTrial2Plot == NaN

if strcmpi(par.PSTH,'on')
 posE=get(ax(NumParts),'Pos');
 posPSTH=get(ax(end),'Pos');
 set(ax(end),'Pos',[posPSTH(1) posPSTH(2)-0.04 posPSTH(3) posE(2)-posPSTH(2)+0.01]);
end

DXC=par.RasterIDbarDX*2;
if par.RasterIDbarWidth~=0
    if (par.RasterPlotNumtrials)
        set(ax(1:NumParts),'YAxisLocation','right');
    else    
        set(ax(1:NumParts),'YColor','w');
    end
    set(yl,'String','');
    for iP=1:NumParts
        pos = get(ax(iP),'Position');
        axb(iP)=axes('Position',[pos(1)-par.RasterIDbarDX-DXC pos(2) DXC pos(4)]);
        %plot([0  0],[0 1],'LineWidth',par.RasterIDbarWidth,'LineColor',par.Colors{iP});
        jP = NewIndex(iP);
        pb(iP)=patch([0 0 0.5 0.5],[1 0 0 1],par.Colors{jP});
        if isfield(par,'Colors2')
            set(pb(iP),'EdgeColor',par.Colors2{jP},'LineWidth',2);
        else
            set(pb(iP),'EdgeColor','none');
        end
        axis fill
        axis off
    end
end

function child_ax=splitaxes(fhandle0,numchilds)
% splits an axis into multiple axes.
parent_ax=get(fhandle0,'Position');
width=parent_ax(:,4)/numchilds;
newaxpos=repmat(parent_ax,numchilds,1);
newaxpos(:,4)=width;
newaxpos(:,2)=fliplr(newaxpos(1,2)+[0:length(newaxpos(:,2))-1]*width);
% child_ax=zeros(1,numchilds);
for iA=1:numchilds,
    child_ax(iA)=subplot('Position',newaxpos(iA,:));
end