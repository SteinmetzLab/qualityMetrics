function qcFigSave(figType, figName)
%QCFIGSAVE is a helper function to save out QC pipeline figures.
%
% INPUT
%   figType (char | string)      = type of QC fig (e.g. AP, LFP, etc.)
%   figName (char | string)      = name of QC plot
%
% OUTPUT 
%   n/a

    p = inputParser;
    p.addRequired('figType', @(s)ischar(s)||isstring(s));
    p.addRequired('figName', @(s)ischar(s)||isstring(s));
    p.parse(figType, figName);
    o = p.Results;

    % create save directory in working directory
    saveDir = fullfile(pwd,'figures',o.figType);

    if ~isfolder(saveDir)
        mkdir(saveDir); 
    end

    % save figure
    fname = fullfile(saveDir, o.figName);
    print(fname, '-dpng', '-r150');
    fprintf('saved %s\n', fname);
end