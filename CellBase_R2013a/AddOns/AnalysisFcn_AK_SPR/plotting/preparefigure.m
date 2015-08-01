function fhandles=preparefigure(NumPlots,figuretype,varargin)
% PREPAREFIGURE makes a figure layout and returns handles

switch figuretype
    
    case '1rnp'
        nrows=NumPlots+1;
        ncols=1;
        nplots=nrows*ncols;
        figpos=[552   220   610   682];
        figh=figure(1);
        wm=0.1;
        hm=0.05;
        ipm=0.04;
        rasht=0.3;
        
        psth_ht=((1-(rasht+hm)-(nplots-1)*ipm)-hm)/(nplots-1);
        
        set(gcf,'Position',figpos);
        
        fhandles(1)=subplot('Position',[wm 1-rasht-hm 1-2*wm rasht]);
        co=1;
        for iH=2:nplots,
            bot=(1-rasht-hm)-co*(psth_ht+ipm);co=co+1;
            fhandles(iH)=subplot('Position',[wm bot 1-2*wm psth_ht]);
        end
        
    case '4r4p'
        nrows=2;
        ncols=2;
        nplots=nrows*ncols;
        figpos=[552   220   610   682];
        figh=figure(1);
        set(figh,'Position',figpos)
        wm=0.1;
        hm=0.05;
        ipm=0.02;
        rasht=0.3;
        psth_ht=((1-(rasht+hm)-(nplots-1)*ipm)-hm)/(nplots-1);
        
    case 'nrnp'
        nplots=NumPlots+1;
        figpos=[359    -2   644   686];
        figh=figure(1);
        set(figh,'Position',figpos,0.1,0.075)
        fhandles=set_subplots(nplots,1);
        
    case 'nrnpw'
        nrows=2;
        ncols=2;
        nplots=nrows*ncols;
        figpos=[359    -2   644   686];
        figh=figure(1);
        set(figh,'Position',figpos)
        for iH=1:nplots,
            fhandles(iH)=subplot(nrows,ncols,iH);
        end
end