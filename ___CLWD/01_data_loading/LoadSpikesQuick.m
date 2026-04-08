function [S, cell_names, params] = LoadSpikesQuick(varargin)
% LoadSpikesQuick  Lightweight spike loader for current session directory
%
%   Loads all spike trains (.t files) from the current directory using
%   LoadSpikes. Subtracts session start time from all spike timestamps
%   so that spikes are in session-relative time (matching LoadSessionData
%   convention and SortBrainstemEventLabels3 laser times).
%
%   Returns the spike train struct and a cell array of clean cell name
%   strings (e.g., 'M222-2022-01-21-TT01-1').
%
%   Does NOT load CSCs, eye tracking, saccade data, or any other session
%   data. Use this when you only need spikes.
%
%   INPUTS:
%       varargin    - name-value pairs:
%           'load_questionable_cells', 0/1  (default: 0)
%
%   OUTPUTS:
%       S           - spike train struct (ts object), times in session-
%                     relative seconds (start_time subtracted)
%       cell_names  - cell array of clean cell name strings
%       params      - parameter struct (includes .starttime)
%
%   USAGE:
%       cd('C:\Jeff\U01\datatouse\M222\M222-2022-01-21')
%       [S, cell_names] = LoadSpikesQuick();
%
%       % Select a specific cell:
%       myCell = SelectTS([], S, 3);
%
%   DEPENDENCIES:
%       LoadSpikes, LoadEvents, process_varargin
%

%% Defaults
load_questionable_cells = 0;

%% Process varargin
process_varargin(varargin);

%% Get session start time from Events file
events_ts = LoadEvents([]);
start_label_idx = [];
for iL = 1:length(events_ts.label)
    if strcmp(events_ts.label{iL}, 'Starting Recording')
        start_label_idx = iL;
        break
    end
end
if isempty(start_label_idx)
    error('Could not find ''Starting Recording'' event in Events file.');
end
starttime = events_ts.t{start_label_idx}(1);  % first start timestamp

%% Load spikes
cfg_spk = [];
cfg_spk.load_questionable_cells = load_questionable_cells;
S = LoadSpikes(cfg_spk);

%% Subtract start time (convert to session-relative seconds)
for iC = 1:length(S.t)
    S.t{iC} = S.t{iC} - starttime;
end

%% Build clean cell name list
nCells = length(S.label);
cell_names = cell(1, nCells);
for iC = 1:nCells
    [~, cell_names{iC}] = fileparts(S.label{iC});
end

%% Store params
params = struct();
params.function_name            = 'LoadSpikesQuick';
params.version                  = '1.1';
params.date_run                 = datestr(now, 'yyyy-mm-dd');
params.load_questionable_cells  = load_questionable_cells;
params.starttime                = starttime;

end

%% Version History
% 2026-03-18  v1.1
%                   - Subtracts session start time from spike timestamps
%                     to match LoadSessionData / SortBrainstemEventLabels3
%                     convention (session-relative seconds).
% 2026-03-18  v1.0 Initial version.
