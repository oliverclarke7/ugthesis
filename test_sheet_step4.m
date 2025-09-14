%% export_EPF1_metrics_singleRowPerRat.m
% This script computes the explore latency for each session.
% Explore latency is defined as the difference between the candidate CP trial and the reversal trial.
% A candidate CP is the first decrease in CP slope (relative to the previous CP) within a window starting at
% (reversal_trial - maxTrialsBeforeReversal) and extending past the reversal.
% Optionally, if exploreSlopeThreshold is specified, the candidate CP's slope must be below that threshold.
%
% After computing session-level explore latency, the script averages the metric within each rat (grouped by rat and age).
% Finally, it pivots the results so that each rat appears on one row and outputs a spreadsheet with columns for each age:
% e.g., EPF1_P30, EPF1_P50, etc.
%
% Adjustable parameters are defined below.

clear; clc;

%% PARAMETERS
critVal = 3.2;                  % Critical value for file selection
maxTrialsBeforeReversal = 10;   % Maximum number of trials BEFORE reversal to consider CPs
exploreSlopeThreshold = [];     % (Optional) Slope threshold; set to [] or Inf to ignore
groupsToInclude = {'SOC'};      % Specify which group(s) to include (set [] to include all)

% File locations
dataDir = '/Users/olly/Desktop/ASI/ASI/MAT Files/';
dataFile = sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal);
filePath = fullfile(dataDir, dataFile);
outputFile = fullfile(dataDir, 'EPF1_test.csv');

%% Load Data
data = readtable(filePath);
data.Properties.VariableNames = lower(data.Properties.VariableNames);

% If groupsToInclude is provided, filter the data accordingly.
if ~isempty(groupsToInclude)
    data = data(ismember(data.group, groupsToInclude), :);
end

% Get unique ages and groups (for later use)
ages = unique(data.age);
groups = unique(data.group);

%% Compute Session-Level Explore Latency
nSessions = height(data);
exploreLatency = NaN(nSessions, 1);

for i = 1:nSessions
    reversalTrial = data.reversal_trial(i);
    
    % Get CP trials and slopes (both before and after reversal)
    cp_before_trials = [data.cp1_trial_before(i), data.cp2_trial_before(i), data.cp3_trial_before(i), data.cp4_trial_before(i)];
    cp_before_slopes = [data.cp1_slope_before(i), data.cp2_slope_before(i), data.cp3_slope_before(i), data.cp4_slope_before(i)];
    cp_after_trials = [data.cp1_trial_after(i), data.cp2_trial_after(i), data.cp3_trial_after(i), data.cp4_trial_after(i)];
    cp_after_slopes = [data.cp1_slope_after(i), data.cp2_slope_after(i), data.cp3_slope_after(i), data.cp4_slope_after(i)];
    
    % Combine CPs into one list
    cp_all_trials = [cp_before_trials, cp_after_trials];
    cp_all_slopes = [cp_before_slopes, cp_after_slopes];
    
    % Remove entries with NaN trial or slope
    validIdx = ~isnan(cp_all_trials) & ~isnan(cp_all_slopes);
    cp_all_trials = cp_all_trials(validIdx);
    cp_all_slopes = cp_all_slopes(validIdx);
    
    % Sort CPs by trial number (ascending)
    [sortedTrials, sortIdx] = sort(cp_all_trials);
    sortedSlopes = cp_all_slopes(sortIdx);
    
    % Define window: consider CPs at or after (reversal - maxTrialsBeforeReversal)
    window_start = reversalTrial - maxTrialsBeforeReversal;
    windowIdx = sortedTrials >= window_start;
    sortedTrials = sortedTrials(windowIdx);
    sortedSlopes = sortedSlopes(windowIdx);
    
    % Initialize candidate latency as NaN (if no candidate is found)
    candidateLatency = NaN;
    
    % Look for the first decrease in slope (starting from the second CP in the window)
    for j = 2:length(sortedTrials)
        if sortedSlopes(j) < sortedSlopes(j-1)
            % If a threshold is provided, check that the candidate CP's slope is below it.
            if ~isempty(exploreSlopeThreshold) && ~isinf(exploreSlopeThreshold)
                if sortedSlopes(j) >= exploreSlopeThreshold
                    continue;  % Not low enough; skip to next candidate.
                end
            end
            % Compute latency relative to reversal (negative means before reversal)
            candidateLatency = sortedTrials(j) - reversalTrial;
            break;
        end
    end
    
    exploreLatency(i) = candidateLatency;
end

% Add the computed explore latency to the data table
data.explore_latency = exploreLatency;

%% Within-Rat Averaging: Compute Rat-Level Explore Latency
% (Assumes a "rat" column exists for subject identification.)
[G, ratID, ratAge, ratGroup] = findgroups(data.rat, data.age, data.group);
ratExploreLatency = splitapply(@(x) mean(x, 'omitnan'), data.explore_latency, G);

% Build a table of rat-level averages (each row corresponds to a unique (rat, age))
T_rat = table(ratID, ratAge, ratGroup, ratExploreLatency);

%% Pivot: One Row Per Rat with Columns for Each Age
uniqueRats = unique(T_rat.ratID, 'stable');  % Maintain order
uniqueAges = unique(T_rat.ratAge);
nRats = numel(uniqueRats);
nAges = numel(uniqueAges);

% Initialize matrix for the explore latency metric
EPF1 = NaN(nRats, nAges);

% Fill the matrix: for each rat, place the explore latency for each age in the proper column.
for iRat = 1:nRats
    currRat = uniqueRats{iRat};  % currRat is a string
    % Use strcmp to find rows for the current rat
    idxRat = strcmp(T_rat.ratID, currRat);
    rowsForRat = T_rat(idxRat, :);
    for k = 1:height(rowsForRat)
        thisAge = rowsForRat.ratAge(k);
        ageCol = find(uniqueAges == thisAge);
        EPF1(iRat, ageCol) = rowsForRat.ratExploreLatency(k);
    end
end

% Create variable names for each age column (e.g., EPF1_P30)
colNames = arrayfun(@(a) sprintf('EPF1_P%d', a), uniqueAges, 'UniformOutput', false);

% Convert the matrix to a table and add the rat IDs
T_EPF1 = array2table(EPF1, 'VariableNames', colNames);
T_out = table(uniqueRats, 'VariableNames', {'rat'});
T_out = [T_out, T_EPF1];

%% Write the output table to CSV
writetable(T_out, outputFile);
fprintf('Exported %d rows (one per rat) to %s\n', height(T_out), outputFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local helper function for computing mean ignoring NaNs
function y = localNanMean(x, dim)
    if nargin < 2
        dim = find(size(x)~=1, 1);
        if isempty(dim)
            dim = 1;
        end
    end
    sumX = sum(x, dim, 'omitnan');
    countX = sum(~isnan(x), dim);
    y = sumX ./ countX;
    y(countX == 0) = NaN;
end