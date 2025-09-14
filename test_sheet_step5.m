%% export_EPF2_metrics_singleRowPerRat.m
% This script computes three explore CP slope metrics for each session:
%   1. Absolute decrease in slope
%   2. Percentage decrease in slope
%   3. Candidate slope (the slope after the decrease)
%
% The explore CP is defined as the first CP (from a combined list of CPs before and after reversal,
% restricted to those occurring at or after (reversal_trial - maxTrialsBeforeReversal))
% that shows a decrease relative to the immediately preceding CP.
%
% After computing session-level metrics, the script averages them within each rat (grouped by rat and age),
% then pivots the results so that each rat appears on one row with columns for each metric by age.
% The output CSV has columns named:
%   EPF2_1_P<age> for absolute decrease,
%   EPF2_2_P<age> for percentage decrease,
%   EPF2_3_P<age> for candidate slope.
%
% Adjustable parameters are defined below.

clear; clc;

%% PARAMETERS
critVal = 3.2;                   % Critical value for file selection
maxTrialsBeforeReversal = 10;    % Maximum trials BEFORE reversal to consider CPs
exploreSlopeThreshold = [];      % (Optional) Slope threshold; set to [] or Inf to ignore
groupsToInclude = {'SOC'};       % Specify which group(s) to include. Set [] to include all.

% File locations
dataDir = '/Users/olly/Desktop/ASI/ASI/MAT Files/';
dataFile = sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal);
filePath = fullfile(dataDir, dataFile);
outputFile = fullfile(dataDir, 'EPF2_test.csv');

%% Load Data
data = readtable(filePath);
data.Properties.VariableNames = lower(data.Properties.VariableNames);

% Filter by group if specified
if ~isempty(groupsToInclude)
    data = data(ismember(data.group, groupsToInclude), :);
end

% Get unique ages and groups (for later pivoting/averaging)
ages = unique(data.age);
groups = unique(data.group);

nSessions = height(data);

%% Compute Session-Level Explore CP Slope Metrics

% Preallocate vectors
absDecrease = NaN(nSessions, 1);     % Absolute decrease in slope
pctDecrease = NaN(nSessions, 1);       % Percentage decrease in slope
candidateSlope = NaN(nSessions, 1);    % Candidate slope (what it decreased to)

for i = 1:nSessions
    reversalTrial = data.reversal_trial(i);
    
    % Get CPs before and after reversal
    cp_before_trials = [data.cp1_trial_before(i), data.cp2_trial_before(i), ...
                        data.cp3_trial_before(i), data.cp4_trial_before(i)];
    cp_before_slopes = [data.cp1_slope_before(i), data.cp2_slope_before(i), ...
                        data.cp3_slope_before(i), data.cp4_slope_before(i)];
    cp_after_trials = [data.cp1_trial_after(i), data.cp2_trial_after(i), ...
                       data.cp3_trial_after(i), data.cp4_trial_after(i)];
    cp_after_slopes = [data.cp1_slope_after(i), data.cp2_slope_after(i), ...
                       data.cp3_slope_after(i), data.cp4_slope_after(i)];
    
    % Combine CPs into one list
    cp_all_trials = [cp_before_trials, cp_after_trials];
    cp_all_slopes = [cp_before_slopes, cp_after_slopes];
    
    % Remove any entries with NaN in trial or slope
    validIdx = ~isnan(cp_all_trials) & ~isnan(cp_all_slopes);
    cp_all_trials = cp_all_trials(validIdx);
    cp_all_slopes = cp_all_slopes(validIdx);
    
    % Sort CPs in ascending order by trial number
    [sortedTrials, sortIdx] = sort(cp_all_trials);
    sortedSlopes = cp_all_slopes(sortIdx);
    
    % Restrict to CPs occurring at or after (reversal_trial - maxTrialsBeforeReversal)
    window_start = reversalTrial - maxTrialsBeforeReversal;
    windowIdx = sortedTrials >= window_start;
    sortedTrials = sortedTrials(windowIdx);
    sortedSlopes = sortedSlopes(windowIdx);
    
    % Initialize candidate metrics as NaN (if no candidate found)
    candAbs = NaN;
    candPct = NaN;
    candSlope = NaN;
    
    % Look for the first decrease in slope (starting at the second CP in the window)
    for j = 2:length(sortedTrials)
        if sortedSlopes(j) < sortedSlopes(j-1)
            % If a threshold is provided, check that the candidate's slope is below it.
            if ~isempty(exploreSlopeThreshold) && ~isinf(exploreSlopeThreshold)
                if sortedSlopes(j) >= exploreSlopeThreshold
                    continue;  % Does not meet threshold, so skip.
                end
            end
            % Compute metrics:
            candAbs = sortedSlopes(j-1) - sortedSlopes(j);    % Absolute decrease
            if sortedSlopes(j-1) > 0
                candPct = (candAbs / sortedSlopes(j-1)) * 100;  % Percentage decrease
            else
                candPct = NaN;
            end
            candSlope = sortedSlopes(j);  % Candidate slope (what it decreased to)
            break;  % Only consider the first decrease
        end
    end
    
    absDecrease(i) = candAbs;
    pctDecrease(i) = candPct;
    candidateSlope(i) = candSlope;
end

% Append session-level metrics to the data table
data.abs_decrease = absDecrease;
data.pct_decrease = pctDecrease;
data.candidate_slope = candidateSlope;

%% Within-Rat Averaging: Compute Rat-Level Averages
% (Assumes a "rat" column exists for subject identification)
[G, ratID, ratAge, ratGroup] = findgroups(data.rat, data.age, data.group);
ratAbsDecrease = splitapply(@(x) mean(x, 'omitnan'), data.abs_decrease, G);
ratPctDecrease = splitapply(@(x) mean(x, 'omitnan'), data.pct_decrease, G);
ratCandidateSlope = splitapply(@(x) mean(x, 'omitnan'), data.candidate_slope, G);

% Build a table of rat-level averages (each row corresponds to a unique (rat, age))
T_rat = table(ratID, ratAge, ratGroup, ratAbsDecrease, ratPctDecrease, ratCandidateSlope);

%% Pivot Data: One Row Per Rat, Columns for Each Metric by Age
uniqueRats = unique(T_rat.ratID, 'stable');  % Keep original ordering
uniqueAges = unique(T_rat.ratAge);
nRats = numel(uniqueRats);
nAges = numel(uniqueAges);

% Preallocate matrices for each metric
EPF2_1 = NaN(nRats, nAges);  % For absolute decrease (metric 1)
EPF2_2 = NaN(nRats, nAges);  % For percentage decrease (metric 2)
EPF2_3 = NaN(nRats, nAges);  % For candidate slope (metric 3)

% Fill the matrices: for each rat, assign the value for each age
for iRat = 1:nRats
    currRat = uniqueRats{iRat};  % rat id is a string
    % Use strcmp for string matching
    idxRat = strcmp(T_rat.ratID, currRat);
    rowsForRat = T_rat(idxRat, :);
    for k = 1:height(rowsForRat)
        thisAge = rowsForRat.ratAge(k);
        ageCol = find(uniqueAges == thisAge);
        EPF2_1(iRat, ageCol) = rowsForRat.ratAbsDecrease(k);
        EPF2_2(iRat, ageCol) = rowsForRat.ratPctDecrease(k);
        EPF2_3(iRat, ageCol) = rowsForRat.ratCandidateSlope(k);
    end
end

% Create variable names for each age column:
% For absolute decrease: EPF2_1_P<age>
% For percentage decrease: EPF2_2_P<age>
% For candidate slope: EPF2_3_P<age>
colNames_1 = arrayfun(@(a) sprintf('EPF2_1_P%d', a), uniqueAges, 'UniformOutput', false);
colNames_2 = arrayfun(@(a) sprintf('EPF2_2_P%d', a), uniqueAges, 'UniformOutput', false);
colNames_3 = arrayfun(@(a) sprintf('EPF2_3_P%d', a), uniqueAges, 'UniformOutput', false);

% Convert matrices to tables
T_EPF2_1 = array2table(EPF2_1, 'VariableNames', colNames_1);
T_EPF2_2 = array2table(EPF2_2, 'VariableNames', colNames_2);
T_EPF2_3 = array2table(EPF2_3, 'VariableNames', colNames_3);

% Build the final output table with one row per rat
T_out = table(uniqueRats, 'VariableNames', {'rat'});
T_out = [T_out, T_EPF2_1, T_EPF2_2, T_EPF2_3];

%% Write the Output Table to CSV
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