function [ccw_onsets, cw_onsets, chenhao_struct] = load_chenhao_saccades(SSN, varargin)
% load_chenhao_saccades  Load Chenhao's saccade onset times for a session.
%
%   Loads saccade_ONSET_all_unit.mat and extracts CCW (temporal) and CW
%   (nasal) saccade onset times for the specified session. NaN values
%   are removed automatically.
%
%   For batch use, pass the loaded struct back in via 'chenhao_struct'
%   to avoid reloading the file on every call.
%
%   INPUTS:
%       SSN             - session name (e.g., 'M039-2020-08-21-1')
%                         Pass '' to load the file without extracting
%                         (returns empty onset arrays, populated struct).
%       varargin        - name-value pairs:
%           chenhao_sacc_path  - path to .mat file
%                                (default: C:\Jeff\U01\datatouse\saccade_ONSET_all_unit.mat)
%           chenhao_struct     - preloaded struct (skips file load)
%
%   OUTPUTS:
%       ccw_onsets      - column vector of CCW (temporal) onset times (s)
%       cw_onsets       - column vector of CW (nasal) onset times (s)
%       chenhao_struct  - the loaded struct (pass back in for batch use)
%
%   USAGE:
%       % Single neuron:
%       [ccw, cw] = load_chenhao_saccades('M079-2020-10-28');
%
%       % Batch mode — load once, reuse struct:
%       [~, ~, cs] = load_chenhao_saccades('');
%       [ccw1, cw1] = load_chenhao_saccades('M079-2020-10-28', 'chenhao_struct', cs);
%       [ccw2, cw2] = load_chenhao_saccades('M222-2022-01-21', 'chenhao_struct', cs);
%
%   NOTES:
%       - temporalSaccades_ONSET = CCW (temporal = away from nose)
%       - nasalSaccades_ONSET    = CW  (nasal = toward nose)
%       - Times are session-relative (same coordinate system as spike times
%         from LoadSpikesQuick)
%       - Indexed by session (1:54 in NPH_list), not by neuron (1:104)
%
%   DEPENDENCIES:
%       process_varargin

%% Defaults
chenhao_sacc_path = 'C:\Jeff\U01\datatouse\saccade_ONSET_all_unit.mat';
chenhao_struct    = [];

%% Process varargin
process_varargin(varargin);

%% Load file if not preloaded
if isempty(chenhao_struct)
    chenhao_struct = load(chenhao_sacc_path, 'NPH_list', ...
        'nasalSaccades_ONSET', 'temporalSaccades_ONSET');
end

%% Handle load-only mode (empty SSN)
if isempty(SSN)
    ccw_onsets = [];
    cw_onsets  = [];
    return
end

%% Find session index
sess_idx = [];
for ii = 1:length(chenhao_struct.NPH_list)
    if strcmp(chenhao_struct.NPH_list{ii}, SSN)
        sess_idx = ii;
        break
    end
end

if isempty(sess_idx)
    error('load_chenhao_saccades:sessionNotFound', ...
        'Session ''%s'' not found in Chenhao saccade onset file (%d sessions available).', ...
        SSN, length(chenhao_struct.NPH_list));
end

%% Extract and clean
ccw_onsets = chenhao_struct.temporalSaccades_ONSET{sess_idx}(:);  % temporal = CCW
cw_onsets  = chenhao_struct.nasalSaccades_ONSET{sess_idx}(:);     % nasal = CW
ccw_onsets = ccw_onsets(~isnan(ccw_onsets));
cw_onsets  = cw_onsets(~isnan(cw_onsets));

end

%% Version History
% 2026-04-05  v1.0  Initial version. Extracted from
%                    analyze_stationary_saccade_bursts_ver_1_4 and
%                    plot_stationary_bout_snippets_ver_1_1 to eliminate
%                    duplicated Chenhao saccade loading logic.
