%%%%%%%%%%%%%%%%%
%%  CellBase   %%
%%   ver 0.5   %%
%%  AK 11/2006 %%
%%%%%%%%%%%%%%%%%
%
%
% Startup functions
%
%   initcb          -   initializes path & builds database (needs to be set once)
%   prealignSpikes  -   sorts spikes into trials and references to events (run once)
%   check           -   checks if everything is consistent
%   loadcb          -   loads cellbase files (main, spikes or events)
% 
% Management functions
%
%   addanalysis   -   adds an analysis to all cells 
%   delanalysis   -   deletes an analysis from cellbase 
%   findanalysis  -   locates analysis already exists 
%   runanalysis   -   runs an analysis function without adding
%
%   addcell       -   adds a cell to cellbase and performs all analyses 
%   delcell       -   deletes a cell from cellbase 
%   findcell      -   locates a cell in cellbase 
% 
%   insertdata    -   inserts data without an analysis function. 
%   exportdata    -   exports hand entered data from CellBase
%
%   addnewsession   - converts Events and prealigns spikes for new sessions
%   addnewcells     - finds & adds new cells and runs all analyses
%
% Search functions
%
%   listtag       - lists cellbase tags 
%   listcell      - lists all properties of a cell  
%   getvalue      - gets values associated with a property for each cell 
%
%   selectcell    -  Returns the cellids satisfying the required criteria
%   selecttrial   -  Returns the trial numbers satisfying the required criteria
%
% Analysis functions
%
%   meanrate                  - example analysis function  
%   calc_selectivity_epoch    - rocAREA selectivity index comparing two epcohs
%   calc_selectivity_side     - rocAREA selectivity index for sides in a given epoch
%   calc_selectivity_outcome  - rocAREA selectivity index for outcome in a given epoch
%
%
% ---------------------------------------------------------------------
%
% To get started check out example_analysis.
%
% For even more fun type 'help cellbase2'
