% start_of_session_copy.m
%
% Run this at the start of every Claude session to populate ___CLWD\for_claude
% with all files Claude needs: Claude-written code + all mvdmlab dependencies.
%
% After running, upload the entire contents of for_claude\ to Claude.
%
% Last verified: 2026-03-17 (all 30 files OK)
% =========================================================================

base     = 'C:\Jeff\Github_repositories\HD-system-brainstem\';
vdm_base = 'C:\Jeff\Github_repositories\vandermeerlab\code-matlab\shared\';
dst      = fullfile(base, '___CLWD', 'for_claude');

if exist(dst, 'dir')
    delete(fullfile(dst, '*'));
else
    mkdir(dst);
end

files_to_copy = {
    % ── Claude-written (___CLWD\) ──────────────────────────────────────
    fullfile(base, '___CLWD',           'plot_AHV_vs_EV_figure.m');
    fullfile(base, '___CLWD',           'plot_EV_TC_from_workspace.m');
    fullfile(base, '___CLWD',           'screen_AHV_vs_EV_opposition.m');
    fullfile(base, '___CLWD',           'diagnose_desaccade_order.m');
    % ── NPH slow phase\ ───────────────────────────────────────────────
    fullfile(base, 'slow phase',         'getSlowPhaseData2.m');
    fullfile(base, 'slow phase',         'getSlowPhaseTC3.m');
    % ── NPH helper functions\ ─────────────────────────────────────────
    fullfile(base, 'helper functions',   'process_varargin.m');
    fullfile(base, 'helper functions',   'HD_GetSSN.m');
    fullfile(base, 'helper functions',   'LoadSpikesJeff.m');
    fullfile(base, 'helper functions',   'LoadSpikesJeff2.m');
    % ── NPH __in process\ ────────────────────────────────────────────
    fullfile(base, '__in process',       'LoadSessionData.m');
    % ── vandermeerlab proc\ ───────────────────────────────────────────
    fullfile(vdm_base, 'proc',           'TuningCurves.m');
    fullfile(vdm_base, 'proc',           'MakeQfromS.m');
    % ── vandermeerlab util\ ───────────────────────────────────────────
    fullfile(vdm_base, 'util',           'nearest_idx3.m');
    fullfile(vdm_base, 'util',           'restrict.m');
    fullfile(vdm_base, 'util',           'restrict2.m');
    fullfile(vdm_base, 'util',           'restrict_idx.m');
    fullfile(vdm_base, 'util',           'antirestrict.m');
    fullfile(vdm_base, 'util',           'trim_histc.m');
    fullfile(vdm_base, 'util',           'histcn.m');
    % ── vandermeerlab datatypes\ ──────────────────────────────────────
    fullfile(vdm_base, 'datatypes\tsd',  'tsd.m');
    fullfile(vdm_base, 'datatypes\tsd',  'CheckTSD.m');
    fullfile(vdm_base, 'datatypes\ts',   'SelectTS.m');
    fullfile(vdm_base, 'datatypes\ts',   'CheckTS.m');
    % ── vandermeerlab io\ ─────────────────────────────────────────────
    fullfile(vdm_base, 'io',             'ProcessConfig.m');
    fullfile(vdm_base, 'io',             'ProcessConfig2.m');
    fullfile(vdm_base, 'io',             'LoadExpKeys.m');
    fullfile(vdm_base, 'io',             'LoadSpikes.m');
    fullfile(vdm_base, 'io',             'pushdir.m');
    fullfile(vdm_base, 'io',             'popdir.m');
};

n_ok   = 0;
n_miss = 0;
for iF = 1:numel(files_to_copy)
    if exist(files_to_copy{iF}, 'file')
        copyfile(files_to_copy{iF}, dst);
        fprintf('OK:        %s\n', files_to_copy{iF});
        n_ok = n_ok + 1;
    else
        fprintf('NOT FOUND: %s\n', files_to_copy{iF});
        n_miss = n_miss + 1;
    end
end

fprintf('\n%d OK, %d not found.\n', n_ok, n_miss);
fprintf('Done — upload contents of for_claude\\ to Claude.\n');
