function switchcb
%SWITCHCB   Switch CellBase data path and file name. 
%   SWITCHCB resets CellBase datapath and file name. For switching between
%   already initialized instances of CellBase, you can also use CHOOSECB.
%
%   See also CHOOSECB.

%   Edit log: BH 5/30/11

% Check current status of CellBase
if isempty(getpref('cellbase'))
    disp('CellBase not initialized. Use initcb first.')
    return
else     % we have some info stored, let's find out how much
    fields = fieldnames(getpref('cellbase'));
    if ismember('datapath',fields)
        datapath = getpref('cellbase','datapath');
    else
        datapath = default_path;
    end
    if ismember('fname',fields)
        fname = getpref('cellbase','fname');
    else
        fname = default_fname;
    end
end

% Reset CellBase directory and name
d = uigetdir(datapath,'Select the root data directory');
if d ~= 0
    datapath = d;
    setpref('cellbase','datapath',datapath);
    %this may be annoying functionality
    %changes the path of fname to the new path
    [pathstr,name,ext,v] = fileparts(fname);
    fname = fullfile(datapath,[name ext]);
    f = uigetfile({'*.mat','MAT-files (*.mat)'},'Select the main database file',fname);
    if f ~=0
        fname = fullfile(datapath,f);
        if exist(fname,'file')
            setpref('cellbase','fname',fname);
        else
            disp('SWITCHCB: No such database file.')
            return
        end
    else
        disp('SWITCHCB canceled');
        return
    end
else
    disp('SWITCHCB canceled');
    return
end

disp(['Data path is ' datapath]);
disp(['File name is ' fname])