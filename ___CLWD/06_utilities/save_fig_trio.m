function save_fig_trio(fh, filepath_base, varargin)
% save_fig_trio  Save a figure as .png, .fig, and .eps
%
%   Saves a figure handle to up to three formats. Handles the known
%   MATLAB savefig bug with invisible figures via try-catch that
%   temporarily toggles visibility. Closes the figure after saving
%   unless doShow is true.
%
%   INPUTS:
%       fh              - figure handle
%       filepath_base   - base path without extension
%                         (e.g., 'C:\...\figures\my_plot')
%       varargin        - name-value pairs:
%           doSavePNG   - save .png (default: 1)
%           doSaveFIG   - save .fig (default: 1)
%           doSaveEPS   - save .eps (default: 1)
%           doShow      - if 0, closes figure after saving (default: 0)
%           png_dpi     - resolution for .png (default: 600)
%
%   USAGE:
%       save_fig_trio(fh, fullfile(out_folder, 'my_figure'));
%       save_fig_trio(fh, path, 'doSaveEPS', 0, 'doShow', 1);
%       save_fig_trio(fh, path, 'doShow', doShow);
%
%   DEPENDENCIES:
%       process_varargin

%% Defaults
doSavePNG = 1;
doSaveFIG = 1;
doSaveEPS = 1;
doShow    = 0;
png_dpi   = 600;

%% Process varargin
process_varargin(varargin);

%% Validate
if ~isvalid(fh)
    warning('save_fig_trio:invalidHandle', ...
        'Invalid figure handle — skipping save.');
    return
end

%% Save .png
if doSavePNG
    print(fh, [filepath_base '.png'], '-dpng', sprintf('-r%d', png_dpi));
end

%% Save .fig (with savefig visibility workaround)
if doSaveFIG
    try
        savefig(fh, [filepath_base '.fig']);
    catch
        set(fh, 'Visible', 'on');
        savefig(fh, [filepath_base '.fig']);
        if ~doShow; set(fh, 'Visible', 'off'); end
    end
end

%% Save .eps
if doSaveEPS
    print(fh, [filepath_base '.eps'], '-depsc', '-painters');
end

%% Close if not showing
if ~doShow
    close(fh);
end

end

%% Version History
% 2026-04-05  v1.0  Initial version. Extracted from
%                    analyze_stationary_saccade_bursts_ver_1_4 local
%                    function. Now a shared utility in 06_utilities/.
