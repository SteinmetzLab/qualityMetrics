function outPath = qcSavePDF(section, varargin)
%QCSAVEPDF  Collate a pipeline section's QC figures into one multi-page PDF.
%
%   outPath = qcSavePDF(section, Name, Value, ...)
%
% Writes each open figure (or each handle passed in 'Figs') as one page of a
% single PDF, so a section's QC plots bundle into one reviewable file: one PDF
% for AP, one for LFP, one for post-Kilosort. By default it grabs every open
% figure, so the usual flow is close all, run the section's plots, then call it.
%
% INPUTS
%   section         (char)      = pipeline section. Aliases set a tidy file name:
%                                 'ap'->AP, 'lfp'->LFP, 'postks'->postKS. Any
%                                 other string is used as-is (sanitized).
%
% NAME-VALUE OPTIONS
%   'Figs'          (graphics)  = figure handles to write, in page order
%                                 (default [] = all open figures, by fig number)
%   'SaveDir'       (char)      = output folder     (default fullfile(pwd,'figures'))
%   'FileName'      (char)      = output file name        (default '<Section>_qc.pdf')
%   'Append'        (logical)   = append to an existing PDF of the same name
%                                 (default false = start fresh)
%   'Resolution'    (double)    = dpi for embedded images         (default 150)
%
% OUTPUT
%   outPath         (char)      = full path of the PDF written
%
% See also qcFigSave, exportgraphics.

    %% parameters
    p = inputParser;
    p.addRequired('section', @(s)ischar(s)||isstring(s));
    p.addParameter('Figs',       [],    @(h) isempty(h) || all(isgraphics(h(:),'figure')));
    p.addParameter('SaveDir',    fullfile(pwd,'figures'), @(s)ischar(s)||isstring(s));
    p.addParameter('FileName',   '',    @(s)ischar(s)||isstring(s));
    p.addParameter('Append',     false, @islogical);
    p.addParameter('Resolution', 150,   @isnumeric);
    p.parse(section, varargin{:});
    o = p.Results;

    % ---- resolve section -> tidy label / default file name ----
    label = sectionLabel(char(o.section));
    if isempty(o.FileName)
        fname = sprintf('%s_qc.pdf', label);
    else
        fname = char(o.FileName);
        if ~endsWith(lower(fname), '.pdf'), fname = [fname '.pdf']; end
    end

    % ---- collect figures ----
    if isempty(o.Figs)
        figs = findall(0, 'Type', 'figure');
        if isempty(figs)
            warning('qcSavePDF:noFigures', ...
                'No open figures to save for section "%s".', label);
            outPath = '';
            return;
        end
        % sort by figure number so pages follow creation order
        nums = arrayfun(@figNum, figs);
        [~, ord] = sort(nums);
        figs = figs(ord);
    else
        figs = o.Figs(:);
    end

    % ---- output path ----
    saveDir = char(o.SaveDir);
    if ~isfolder(saveDir), mkdir(saveDir); end
    outPath = fullfile(saveDir, fname);

    % start fresh unless appending
    if ~o.Append && isfile(outPath)
        delete(outPath);
    end

    % append figures
    
    if isempty(which('exportgraphics'))
        error('qcSavePDF:noExportgraphics', ...
            ['exportgraphics not found; multi-page PDF export needs MATLAB ' ...
             'R2020a or later.']);
    end
    nPage = numel(figs);
    for k = 1:nPage
        f = figs(k);
        try
            exportgraphics(f, outPath, 'Append', true, ...
                'ContentType', 'vector', 'Resolution', o.Resolution);
        catch ME
            warning('qcSavePDF:pageFailed', ...
                'Skipped figure %g (page %d/%d): %s', figNum(f), k, nPage, ME.message);
        end
    end

    fprintf('saved %s  (%d page(s))\n', outPath, nPage);
end

%% helpers

function label = sectionLabel(s)
% Map a section string to a tidy file-name label; pass others through sanitized.
    switch lower(strtrim(s))
        case {'ap','apqc','ap_qc'}
            label = 'AP';
        case {'lfp','lf','lfpqc','lfp_qc'}
            label = 'LFP';
        case {'postks','post-ks','ks','post-kilosort','postkilosort','post_ks'}
            label = 'postKS';
        otherwise
            % keep only filename-safe characters
            label = regexprep(s, '[^\w-]', '_');
            if isempty(label), label = 'section'; end
    end
end

function n = figNum(f)
% Figure Number, or a large fallback for numberless (e.g. uifigure) handles.
    n = get(f, 'Number');
    if isempty(n), n = inf; end
end
