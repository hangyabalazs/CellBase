tic;
session_cellids=unique_session_cells;
problem_Clusep_sessions=[];
for iSess=1:length(session_cellids),
    sessid=session_cellids(iSess);
    cd(cellid2fnames(sessid,'Sess'));
    if exist(cellid2fnames(sessid,'wv'))==2, % this has session/cellid conflict
        disp('This session is already done')
        continue
    else
    try
        Create_CQ_File
    catch
        problem_Clusep_sessions=[problem_Clusep_sessions sessid];
    end
    end
end
toc;