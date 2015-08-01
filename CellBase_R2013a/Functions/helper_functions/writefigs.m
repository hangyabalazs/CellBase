function writefigs(H,pdfname)
%WRITEFIGS   Write figures to pdf.
%   WRITEFIGS(H,PDFNAME) appends figure H to the pdf file PDFNAME. H can be
%   either a struct with figure handles, or a single figure handle.
%
%   See also PRINT.
 
% Append to pdf
if isstruct(H)  % H is either a struct with figure handles, or a single fig. handle
    fls = fieldnames(H);
    for fs = 1:length(fls)
        h = H.(fls{fs});  % current handle
        if ishandle(h)
            export_fig(h,'-append',pdfname, '-zbuffer');  % write to pdf
        end
    end
else
    export_fig(H,'-append',pdfname, '-zbuffer');  % write to pdf
end