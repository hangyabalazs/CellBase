function cellbase2
%
% CellBase
%
% Functions used internally by CellBase 
%
% (You don't really need to know this, but might come in handy.) 
%
% findallcells  -   locates all the clustered mat files in the default path 
% findcellpos   -   returns the position of cell in cellbase 
% findprop      -   returns the column containing a property  
%
% fname2cellid  -   creates a cellid from a valid unit data file name
% cellid2fnames -   converts cellid into names of session & spike data file names
% cellid2tags   -   converts cellid into separate tags 
% cellid2vals   -   converts cellid into separate values (basic cell properties) 
%  
% validcellid   -   checks if cellid is valid
% listfiles     -   lists files in directory
% listdir       -   lists directories in directory
% nargvout      -   returns the size of the vector output argument for a function
% unique_cell   -   returns unique cells in array
%
%