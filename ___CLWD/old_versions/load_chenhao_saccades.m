function [temporal_onsets, nasal_onsets, params] = load_chenhao_saccades(SSN, varargin)
% load_chenhao_saccades  Load Chenhao's saccade onset times for a session.
%
%   Chenhao's onset times are aligned to saccade START (not velocity
%   threshold crossing), making them more accurate for PETH alignment.
%
%   INPUTS:
%     SSN             - session name string (e.g. 'M284-2022-05-12-1')
%
%   OUTPUTS:
%     temporal_onsets  - column vector of temporal (CCW) saccade onset times (s)
%     nasal_onsets     - column vector of nasal (CW) saccade onset times (s)
%     params           - struct with metadata
%
%   USAGE:
%     [tON, nON, par] = load_chenhao_saccades('M284-2022-05-12-1');

%% Defaults
mat_path = 'C:\Jeff\U01\datatouse\saccade_ONSET_all_unit.mat';
process_varargin(varargin);

%% Params
params = struct();
params.function_name = 'load_chenhao_saccades';
params.version       = '1.0';
params.mat_path      = mat_path;
params.SSN           = SSN;
params.source        = 'Chenhao Bao (JHU) — saccade onset re-alignment';

%% Load
persistent chenhao_cache
if isempty(chenhao_cache)
    chenhao_cache = load(mat_path, 'NPH_list', 'temporalSaccades_ONSET', 'nasalSaccades_ONSET');
    fprintf('  Loaded Chenhao saccade onsets (%d sessions)\n', length(chenhao_cache.NPH_list));
end

%% Find session index
sess_idx = [];
for iS = 1:length(chenhao_cache.NPH_list)
    if strcmp(chenhao_cache.NPH_list{iS}, SSN)
        sess_idx = iS;
        break
    end
end

if isempty(sess_idx)
    warning('load_chenhao_saccades: SSN ''%s'' not found in saccade_ONSET_all_unit.mat', SSN);
    temporal_onsets = [];
    nasal_onsets    = [];
    params.status   = 'SSN not found';
    return
end

%% Extract
temporal_onsets = chenhao_cache.temporalSaccades_ONSET{sess_idx}(:);
nasal_onsets    = chenhao_cache.nasalSaccades_ONSET{sess_idx}(:);

% Remove NaNs
temporal_onsets = temporal_onsets(~isnan(temporal_onsets));
nasal_onsets    = nasal_onsets(~isnan(nasal_onsets));

params.sess_idx      = sess_idx;
params.n_temporal    = length(temporal_onsets);
params.n_nasal       = length(nasal_onsets);
params.status        = 'OK';

end
