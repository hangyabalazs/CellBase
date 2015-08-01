function nlxcsc2mat2(pname,varargin)
% e.g.: nlxcsc2mat2(pname,'Channels','all') % 'position','Events'
% this replaces the MatlabRead function. It can convert a specified
% continuous channel from Neuralynx to matlab.
% INPUTS:

% path of the data directory
% Continuous channel no. i.e. 1, 2,'all' for all files or 'Events' to convert events file
% if varargin(2) = 'Events' then read Events.nev into mat file
% If nothing is specified it converts all files from chan 4 to 8.
% SPR 2008-07-01
% default_args={...
%     'ifsave', 1;...
%     'change_fname',1;...
%     'Channels',1;...% also takes 'all' and 'Events'
%     'dt',0.001;... % new sampling rate sf=1/dt; 0.001=1000Hz
%     'FieldSelectionArray',[1 0 0 0 1];...% only TimStamps and Samples
%     'ExtractHeader',0;... % No Header
%     'ExtractMode',1;... % All records
%     'ExtractModeArray',[];... % Extract All records
%     'factored_downsampling',0;... % the proper way of decimating if decimating factor is greater than 13
%     };
% SPR 2010-04-08
% rewrite this so that it takes a defineCSCAssignments file to specify the CSC naming conventions etc.

default_args={...
    'ifsave', 1;...
    'change_fname',1;...
    'Channels',0;... % also acceptable is 'all' and 'Events'; '0' is a special case because then nothing corresponds to it and so filename is not changed.
    'dt',0.001;... % new sampling rate sf=1/dt; 0.001=1000Hz
    'FieldSelectionArray',[1 0 0 0 1];...% only TimStamps and Samples
    'ExtractHeader',0;... % No Header
    'ExtractMode',1;... % All records
    'ExtractModeArray',[];... % Extract All records
    'factored_downsampling',0;... % the proper way of decimating if decimating factor is greater than 13
    };
[g, error] = parse_args(default_args,varargin{:});
if ischar(g.Channels),
    switch g.Channels
        case 'Events'
            g.Channels=1000;
        case 'all'
            g.Channels=[1:40 1000 2000];
        case 'position'
            g.Channels=2000;
    end
end

for iC=1:length(g.Channels),
    tic;
    action = sprintf('CSC%s.Ncs',num2str(g.Channels(iC)));
    switch action
        case 'CSC1.Ncs'
            fname = 'CSC1.Ncs';new_fname='CSC0c1';
        case 'CSC2.Ncs'
            fname = 'CSC2.Ncs';new_fname='CSC0c2';
        case 'CSC3.Ncs'
            fname = 'CSC3.Ncs';new_fname='CSC0c3';
        case 'CSC4.Ncs'
            fname = 'CSC4.Ncs';new_fname='CSC0c4';
        case 'CSC5.Ncs'
            fname = 'CSC5.Ncs';new_fname='CSC0c5';
        case 'CSC6.Ncs'
            fname = 'CSC6.Ncs';new_fname='CSC0c6';
        case 'CSC7.Ncs'
            fname = 'CSC7.Ncs';new_fname='CSC0c7';
        case 'CSC8.Ncs'
            fname = 'CSC8.Ncs';new_fname='CSC0c8';
        case 'CSC9.Ncs'
            fname = 'CSC9.Ncs';new_fname='CSC1c1';
        case 'CSC10.Ncs'
            fname = 'CSC10.Ncs';new_fname='CSC1c2';
        case 'CSC11.Ncs'
            fname = 'CSC11.Ncs';new_fname='CSC1c3';
        case 'CSC12.Ncs'
            fname = 'CSC12.Ncs';new_fname='CSC1c4';
        case 'CSC13.Ncs'
            fname = 'CSC13.Ncs';new_fname='CSC2c1';
        case 'CSC14.Ncs'
            fname = 'CSC14.Ncs';new_fname='CSC2c2';
        case 'CSC15.Ncs'
            fname = 'CSC15.Ncs';new_fname='CSC2c3';
        case 'CSC16.Ncs'
            fname = 'CSC16.Ncs';new_fname='CSC2c4';
        case 'CSC17.Ncs'
            fname = 'CSC17.Ncs';new_fname='CSC3c1';
        case 'CSC18.Ncs'
            fname = 'CSC18.Ncs';new_fname='CSC3c2';
        case 'CSC19.Ncs'
            fname = 'CSC19.Ncs';new_fname='CSC3c3';
        case 'CSC20.Ncs'
            fname = 'CSC20.Ncs';new_fname='CSC3c4';
        case 'CSC21.Ncs'
            fname = 'CSC21.Ncs';new_fname='CSC4c1';
        case 'CSC22.Ncs'
            fname = 'CSC22.Ncs';new_fname='CSC4c2';
        case 'CSC23.Ncs'
            fname = 'CSC23.Ncs';new_fname='CSC4c3';
        case 'CSC24.Ncs'
            fname = 'CSC24.Ncs';new_fname='CSC4c4';
        case 'CSC25.Ncs'
            fname = 'CSC25.Ncs';new_fname='CSC5c1';
        case 'CSC26.Ncs'
            fname = 'CSC26.Ncs';new_fname='CSC5c2';
        case 'CSC27.Ncs'
            fname = 'CSC27.Ncs';new_fname='CSC5c3';
        case 'CSC28.Ncs'
            fname = 'CSC28.Ncs';new_fname='CSC5c4';
        case 'CSC29.Ncs'
            fname = 'CSC29.Ncs';new_fname='CSC6c1';
        case 'CSC30.Ncs'
            fname = 'CSC30.Ncs';new_fname='CSC6c2';
        case 'CSC31.Ncs'
            fname = 'CSC31.Ncs';new_fname='CSC6c3';
        case 'CSC32.Ncs'
            fname = 'CSC32.Ncs';new_fname='CSC6c4';
        case 'CSC33.Ncs'
            fname = 'CSC33.Ncs';new_fname='CSC7c1';
        case 'CSC34.Ncs'
            fname = 'CSC34.Ncs';new_fname='CSC7c2';
        case 'CSC35.Ncs'
            fname = 'CSC35.Ncs';new_fname='CSC7c3';
        case 'CSC36.Ncs'
            fname = 'CSC36.Ncs';new_fname='CSC7c4';
        case 'CSC37.Ncs'
            fname = 'CSC37.Ncs';new_fname='CSC7c5';
        case 'CSC38.Ncs'
            fname = 'CSC38.Ncs';new_fname='CSC7c6';
        case 'CSC39.Ncs'
            fname = 'CSC39.Ncs';new_fname='CSC7c7';
        case 'CSC40.Ncs'
            fname = 'CSC40.Ncs';new_fname='CSC7c8';
        case 'CSC1000.Ncs'
            fname = 'Events.nev';new_fname='EVENTS';
            param2 = [1 1 1 1 1];
            param3 = 1;
            param4 = 1;
            param5 = [];
            [Events_TimeStamps, Events_EventIDs, Events_Nttls, Events_Extras, Events_EventStrings, Events_NlxHeader] = Nlx2MatEV([pname fname],param2,param3,param4,param5);
% 
%             [Events_TimeStamps,Events_EventIDs,Events_Nttls,Events_Extras,Events_EventStrings]...
%                 = Nlx2MatEV([pname filesep fname],[1 1 1 1 1],...
%                 g.ExtractMode,g.ExtractModeArray);
            Events_TimeStamps=Events_TimeStamps*1e-6;
            save([pname filesep new_fname],'Events_EventIDs','Events_EventStrings',...
                'Events_Extras','Events_Nttls','Events_TimeStamps');
        case 'CSC2000.Ncs'
            fname = 'VT1.nvt';new_fname='POSITION';
            param2 = [1 1 1 1 1 1];
            param3 = 0;
            param4 = 1;
            param5 = [];
            [TimeStamps, ExtractedX, ExtractedY, ExtractedAngle, Targets, Points] = Nlx2MatVT( [pname filesep fname],param2,param3,param4,param5);
            % [TimeStamps, ExtractedX, ExtractedY, ExtractedAngle, Targets, Points]...
            %= Nlx2MatVT_v4( Filename, FieldSelection, ExtractHeader, ExtractMode, ModeArray );
            TimeStamps=TimeStamps*1e-6;

% save([path '\' fn(1:end-4)],'TS','ScNumbers', 'CellNumbers', 'Params', 'DataPoints','NlxHeader');
    save([pname filesep new_fname],'TimeStamps', 'ExtractedX', 'ExtractedY', 'ExtractedAngle', 'Targets', 'Points');

        otherwise % no match with filename; leave it as is.
            [PATH,fname,EXT,junk] = fileparts(pname);
            new_fname = fname;
            fname = [fname EXT];
            pname = PATH;
            if strmatch(EXT,'.nev','exact'),
                param2 = [1 1 1 1 1];
                param3 = 1;
                param4 = 1;
                param5 = [];
                [Events_TimeStamps, Events_EventIDs, Events_Nttls, Events_Extras, Events_EventStrings, Events_NlxHeader] = Nlx2MatEV([pname filesep fname],param2,param3,param4,param5);
                %
                %             [Events_TimeStamps,Events_EventIDs,Events_Nttls,Events_Extras,Events_EventStrings]...
                %                 = Nlx2MatEV([pname filesep fname],[1 1 1 1 1],...
                %                 g.ExtractMode,g.ExtractModeArray);
                Events_TimeStamps=Events_TimeStamps*1e-6;
                save([pname filesep new_fname],'Events_EventIDs','Events_EventStrings',...
                    'Events_Extras','Events_Nttls','Events_TimeStamps');
                laptime=toc;
                disp(sprintf('%s converted in %d seconds; %d of %d remaining',fname,laptime,iC,length(g.Channels)))
                return
            end
    end
    if iC~=1000,
        try
            [Samples,dt,ts]=convertcsc(pname,fname,new_fname,g);
        catch
        end
    end
    laptime=toc;
    disp(sprintf('%s converted in %d seconds; %d of %d remaining',fname,laptime,iC,length(g.Channels)))
end

function [Samples,dt,ts]=convertcsc(pname,fname,new_fname,g)
[ts,Samples]=Nlx2MatCSC([pname filesep fname],g.FieldSelectionArray,g.ExtractHeader,g.ExtractMode,g.ExtractModeArray); %#ok<NASGU>
dt=(ts(2)-ts(1))*1e-6/size(Samples,1);
ts=ts(1)*1e-6; % seconds
Samples=Samples(:);
dec_factor=floor(g.dt/dt);
dt=dt*dec_factor;
if dec_factor>13,
    if g.factored_downsampling==1,
        dec_factor=factor(dec_factor);
    end
    for iF=1:length(dec_factor),
        Samples=decimate(Samples,dec_factor(iF));
    end
end
%             if g.down
if g.ifsave==1,
    save([pname filesep new_fname],'Samples','ts','dt');
end