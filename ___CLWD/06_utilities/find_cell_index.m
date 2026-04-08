function cell_idx = find_cell_index(spike_cell_names, neuron_name, SSN)
% find_cell_index  Match a neuron name to loaded spike cell names.
%
%   Tries exact match first. If that fails, strips the session prefix
%   from neuron_name and tries a partial match.
%
%   INPUTS:
%       spike_cell_names - cell array of cell name strings (from
%                          LoadSpikesQuick or LoadSpikes)
%       neuron_name      - target neuron (e.g., 'M222-2022-01-21-TT01_1')
%       SSN              - session name (e.g., 'M222-2022-01-21')
%
%   OUTPUTS:
%       cell_idx         - index into spike_cell_names, or [] if not found
%
%   USAGE:
%       [S, cell_names] = LoadSpikesQuick();
%       idx = find_cell_index(cell_names, 'M222-2022-01-21-TT01_1', 'M222-2022-01-21');
%       myCell = SelectTS([], S, idx);
%
%   DEPENDENCIES:
%       none

cell_idx = [];

%% Try exact match
for cc = 1:length(spike_cell_names)
    if strcmp(spike_cell_names{cc}, neuron_name)
        cell_idx = cc;
        return
    end
end

%% Try partial match (strip session prefix)
suffix = strrep(neuron_name, [SSN '-'], '');
for cc = 1:length(spike_cell_names)
    if contains(spike_cell_names{cc}, suffix)
        cell_idx = cc;
        return
    end
end

end

%% Version History
% 2026-04-05  v1.0  Initial version. Extracted from
%                    analyze_stationary_saccade_bursts_ver_1_4 and
%                    plot_stationary_bout_snippets_ver_1_1 to eliminate
%                    duplicated cell matching logic.
