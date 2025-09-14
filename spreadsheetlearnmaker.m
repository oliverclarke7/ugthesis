%% Integrated Regression Export Script
% This script computes three series of measures:
%  1) LF and ETF measures from session-level data.
%  2) Explore measures.
%  3) Reversal averages.
%
% It aggregates session-level values to the rat level, pivots the data
% into wide format, merges with rat metadata, then adds reversal measures.
% The output is exported to a CSV file whose name includes the crit value,
% preserving uppercase Rat IDs (e.g. "DEV005017010").

clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS for LF/ETF Measures
critVal = 2.6;                % Critical value (also used in file names)
slopeThreshold = 0.6;         % Slope threshold for candidate CP (LF/ETF measures)
minTrialsBeforeNextCP = 5;    % Minimum trials gap to the next CP
minTrialsBeforeReversal = 0;  % Minimum trials gap to the reversal
maxTrialsBeforeReversal = [35]; % Maximum gap (can be [] or Inf for no limit)

% Filtering parameters for LF/ETF:
reversalFilter = 'REV1_REV2'; % For LF measures
groupsToInclude = {'SOC'};    % Only include sessions from group 'SOC'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS for Explore Measures
maxTrialsBeforeReversal_EXP = 10;    % Maximum number of trials BEFORE reversal to consider CPs (explore)
exploreSlopeThreshold_EXP = [];      % (Optional) Slope threshold for explore measures; [] or Inf to ignore
groupsToInclude_EXP = {'SOC'};       % Only include sessions from group 'SOC'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data for LF/ETF Measures
dataDir = '/Users/olly/Desktop/ASI/ASI/MAT Files/';
dataFile = sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal);
filePath = fullfile(dataDir, dataFile);

data_all = readtable(filePath);
data_all.Properties.VariableNames = lower(data_all.Properties.VariableNames);

% Filter to include only the specified groups (i.e. only 'SOC')
data_all = data_all(ismember(data_all.group, groupsToInclude), :);

% Convert Rat IDs to uppercase for consistent merging
data_all.rat = upper(strtrim(data_all.rat));

% Split into LF and ETF datasets:
data_LF = data_all(strcmp(data_all.reversal_name, reversalFilter), :);
data_ETF = data_all(~strcmp(data_all.reversal_name, reversalFilter), :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute LF Measures (REV1_REV2 sessions)

% --- LF1: CP_before Latency ---
cpBeforeLatency = NaN(height(data_LF), 1);
for i = 1:height(data_LF)
    slopes_before = [data_LF.cp1_slope_before(i), data_LF.cp2_slope_before(i), ...
                     data_LF.cp3_slope_before(i), data_LF.cp4_slope_before(i)];
    trials_before = [data_LF.cp1_trial_before(i), data_LF.cp2_trial_before(i), ...
                     data_LF.cp3_trial_before(i), data_LF.cp4_trial_before(i)];
    reversalTrial = data_LF.reversal_trial(i);
    candidateIndices = find(slopes_before >= slopeThreshold & ~isnan(slopes_before));
    for j = 1:length(candidateIndices)
        idx = candidateIndices(j);
        candidateTrial = trials_before(idx);
        gapToReversal = reversalTrial - candidateTrial;
        if gapToReversal < minTrialsBeforeReversal, continue; end
        if ~isempty(maxTrialsBeforeReversal) && ~isinf(maxTrialsBeforeReversal)
            if gapToReversal > maxTrialsBeforeReversal, continue; end
        end
        if idx < length(trials_before) && ~isnan(trials_before(idx+1))
            nextEventTrial = trials_before(idx+1);
        else
            nextEventTrial = reversalTrial;
        end
        if (nextEventTrial - candidateTrial) < minTrialsBeforeNextCP, continue; end
        cpBeforeLatency(i) = candidateTrial;
        break;
    end
end
data_LF.cp_before_latency = cpBeforeLatency;

% --- LF2 Measures ---
nSessions = height(data_LF);
rawFirstCandidateSlope = NaN(nSessions, 1);
weightedCandidateSlope = NaN(nSessions, 1);
weightedCandidateSlopeFig3 = NaN(nSessions, 1);

for i = 1:nSessions
    slopes_before = [data_LF.cp1_slope_before(i), data_LF.cp2_slope_before(i), ...
                     data_LF.cp3_slope_before(i), data_LF.cp4_slope_before(i)];
    trials_before = [data_LF.cp1_trial_before(i), data_LF.cp2_trial_before(i), ...
                     data_LF.cp3_trial_before(i), data_LF.cp4_trial_before(i)];
    reversalTrial = data_LF.reversal_trial(i);
    candidateSlopes = [];
    candidateWeights = [];
    firstFound = false;
    for j = 1:length(slopes_before)
        currSlope = slopes_before(j);
        currTrial = trials_before(j);
        if isnan(currSlope) || isnan(currTrial), continue; end
        if currSlope < slopeThreshold, continue; end
        gapToReversal = reversalTrial - currTrial;
        if gapToReversal < minTrialsBeforeReversal, continue; end
        if ~isempty(maxTrialsBeforeReversal) && ~isinf(maxTrialsBeforeReversal)
            if gapToReversal > maxTrialsBeforeReversal, continue; end
        end
        if j < length(trials_before) && ~isnan(trials_before(j+1))
            nextEventTrial = trials_before(j+1);
        else
            nextEventTrial = reversalTrial;
        end
        if (nextEventTrial - currTrial) < minTrialsBeforeNextCP, continue; end
        candidateSlopes(end+1) = currSlope;
        candidateWeights(end+1) = nextEventTrial - currTrial;
        if ~firstFound
            rawFirstCandidateSlope(i) = currSlope;
            firstFound = true;
        end
    end
    if ~isempty(candidateSlopes)
        weightedCandidateSlope(i) = sum(candidateSlopes .* candidateWeights) / sum(candidateWeights);
    end
end
data_LF.raw_first_candidate_slope = rawFirstCandidateSlope;
data_LF.weighted_candidate_slope = weightedCandidateSlope;

for i = 1:nSessions
    slopes_before = [data_LF.cp1_slope_before(i), data_LF.cp2_slope_before(i), ...
                     data_LF.cp3_slope_before(i), data_LF.cp4_slope_before(i)];
    trials_before = [data_LF.cp1_trial_before(i), data_LF.cp2_trial_before(i), ...
                     data_LF.cp3_trial_before(i), data_LF.cp4_trial_before(i)];
    reversalTrial = data_LF.reversal_trial(i);
    firstCandidateIndex = [];
    for j = 1:length(slopes_before)
        currSlope = slopes_before(j);
        currTrial = trials_before(j);
        if isnan(currSlope) || isnan(currTrial), continue; end
        if currSlope < slopeThreshold, continue; end
        gapToReversal = reversalTrial - currTrial;
        if gapToReversal < minTrialsBeforeReversal, continue; end
        if ~isempty(maxTrialsBeforeReversal) && ~isinf(maxTrialsBeforeReversal)
            if gapToReversal > maxTrialsBeforeReversal, continue; end
        end
        if j < length(trials_before) && ~isnan(trials_before(j+1))
            nextEventTrial = trials_before(j+1);
        else
            nextEventTrial = reversalTrial;
        end
        if (nextEventTrial - currTrial) < minTrialsBeforeNextCP, continue; end
        firstCandidateIndex = j;
        break;
    end
    if ~isempty(firstCandidateIndex)
        candidateSlopesFig3 = [];
        candidateWeightsFig3 = [];
        for j = firstCandidateIndex:length(slopes_before)
            if isnan(trials_before(j)), continue; end
            if trials_before(j) >= reversalTrial, continue; end
            slope_j = slopes_before(j);
            if (j < length(trials_before)) && (~isnan(trials_before(j+1))) && (trials_before(j+1) < reversalTrial)
                weight_j = trials_before(j+1) - trials_before(j);
            else
                weight_j = reversalTrial - trials_before(j);
            end
            candidateSlopesFig3(end+1) = slope_j;
            candidateWeightsFig3(end+1) = weight_j;
        end
        if ~isempty(candidateSlopesFig3)
            weightedCandidateSlopeFig3(i) = sum(candidateSlopesFig3 .* candidateWeightsFig3) / sum(candidateWeightsFig3);
        end
    end
end
data_LF.weighted_candidate_slope_fig3 = weightedCandidateSlopeFig3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute ETF Measures (sessions NOT equal to REV1_REV2)

% --- ETF1: CP_before Latency ---
cpBeforeLatency_ETF = NaN(height(data_ETF), 1);
for i = 1:height(data_ETF)
    slopes_before = [data_ETF.cp1_slope_before(i), data_ETF.cp2_slope_before(i), ...
                     data_ETF.cp3_slope_before(i), data_ETF.cp4_slope_before(i)];
    trials_before = [data_ETF.cp1_trial_before(i), data_ETF.cp2_trial_before(i), ...
                     data_ETF.cp3_trial_before(i), data_ETF.cp4_trial_before(i)];
    reversalTrial = data_ETF.reversal_trial(i);
    candidateIndices = find(slopes_before >= slopeThreshold & ~isnan(slopes_before));
    for j = 1:length(candidateIndices)
        idx = candidateIndices(j);
        candidateTrial = trials_before(idx);
        gapToReversal = reversalTrial - candidateTrial;
        if gapToReversal < minTrialsBeforeReversal, continue; end
        if ~isempty(maxTrialsBeforeReversal) && ~isinf(maxTrialsBeforeReversal)
            if gapToReversal > maxTrialsBeforeReversal, continue; end
        end
        if idx < length(trials_before) && ~isnan(trials_before(idx+1))
            nextEventTrial = trials_before(idx+1);
        else
            nextEventTrial = reversalTrial;
        end
        if (nextEventTrial - candidateTrial) < minTrialsBeforeNextCP, continue; end
        cpBeforeLatency_ETF(i) = candidateTrial;
        break;
    end
end
data_ETF.cp_before_latency_ETF = cpBeforeLatency_ETF;

% --- ETF2 Measures ---
nSessions_ETF = height(data_ETF);
rawFirstCandidateSlope_ETF = NaN(nSessions_ETF, 1);
weightedCandidateSlope_ETF = NaN(nSessions_ETF, 1);
weightedCandidateSlopeFig3_ETF = NaN(nSessions_ETF, 1);

for i = 1:nSessions_ETF
    slopes_before = [data_ETF.cp1_slope_before(i), data_ETF.cp2_slope_before(i), ...
                     data_ETF.cp3_slope_before(i), data_ETF.cp4_slope_before(i)];
    trials_before = [data_ETF.cp1_trial_before(i), data_ETF.cp2_trial_before(i), ...
                     data_ETF.cp3_trial_before(i), data_ETF.cp4_trial_before(i)];
    reversalTrial = data_ETF.reversal_trial(i);
    candidateSlopes = [];
    candidateWeights = [];
    firstFound = false;
    for j = 1:length(slopes_before)
        currSlope = slopes_before(j);
        currTrial = trials_before(j);
        if isnan(currSlope) || isnan(currTrial), continue; end
        if currSlope < slopeThreshold, continue; end
        gapToReversal = reversalTrial - currTrial;
        if gapToReversal < minTrialsBeforeReversal, continue; end
        if ~isempty(maxTrialsBeforeReversal) && ~isinf(maxTrialsBeforeReversal)
            if gapToReversal > maxTrialsBeforeReversal, continue; end
        end
        if j < length(trials_before) && ~isnan(trials_before(j+1))
            nextEventTrial = trials_before(j+1);
        else
            nextEventTrial = reversalTrial;
        end
        if (nextEventTrial - currTrial) < minTrialsBeforeNextCP, continue; end
        candidateSlopes(end+1) = currSlope;
        candidateWeights(end+1) = nextEventTrial - currTrial;
        if ~firstFound
            rawFirstCandidateSlope_ETF(i) = currSlope;
            firstFound = true;
        end
    end
    if ~isempty(candidateSlopes)
        weightedCandidateSlope_ETF(i) = sum(candidateSlopes .* candidateWeights) / sum(candidateWeights);
    end
end
data_ETF.raw_first_candidate_slope_ETF = rawFirstCandidateSlope_ETF;
data_ETF.weighted_candidate_slope_ETF = weightedCandidateSlope_ETF;

for i = 1:nSessions_ETF
    slopes_before = [data_ETF.cp1_slope_before(i), data_ETF.cp2_slope_before(i), ...
                     data_ETF.cp3_slope_before(i), data_ETF.cp4_slope_before(i)];
    trials_before = [data_ETF.cp1_trial_before(i), data_ETF.cp2_trial_before(i), ...
                     data_ETF.cp3_trial_before(i), data_ETF.cp4_trial_before(i)];
    reversalTrial = data_ETF.reversal_trial(i);
    firstCandidateIndex = [];
    for j = 1:length(slopes_before)
        currSlope = slopes_before(j);
        currTrial = trials_before(j);
        if isnan(currSlope) || isnan(currTrial), continue; end
        if currSlope < slopeThreshold, continue; end
        gapToReversal = reversalTrial - currTrial;
        if gapToReversal < minTrialsBeforeReversal, continue; end
        if ~isempty(maxTrialsBeforeReversal) && ~isinf(maxTrialsBeforeReversal)
            if gapToReversal > maxTrialsBeforeReversal, continue; end
        end
        if j < length(trials_before) && ~isnan(trials_before(j+1))
            nextEventTrial = trials_before(j+1);
        else
            nextEventTrial = reversalTrial;
        end
        if (nextEventTrial - currTrial) < minTrialsBeforeNextCP, continue; end
        firstCandidateIndex = j;
        break;
    end
    if ~isempty(firstCandidateIndex)
        candidateSlopesFig3 = [];
        candidateWeightsFig3 = [];
        for j = firstCandidateIndex:length(slopes_before)
            if isnan(trials_before(j)), continue; end
            if trials_before(j) >= reversalTrial, continue; end
            slope_j = slopes_before(j);
            if (j < length(trials_before)) && (~isnan(trials_before(j+1))) && (trials_before(j+1) < reversalTrial)
                weight_j = trials_before(j+1) - trials_before(j);
            else
                weight_j = reversalTrial - trials_before(j);
            end
            candidateSlopesFig3(end+1) = slope_j;
            candidateWeightsFig3(end+1) = weight_j;
        end
        if ~isempty(candidateSlopesFig3)
            weightedCandidateSlopeFig3_ETF(i) = sum(candidateSlopesFig3 .* candidateWeightsFig3) / sum(candidateWeightsFig3);
        end
    end
end
data_ETF.weighted_candidate_slope_fig3_ETF = weightedCandidateSlopeFig3_ETF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Explore Measures using the "explore_find_1" Technique
% Load data for explore measures (filter by groupsToInclude_EXP)
data_EXP = readtable(filePath);
data_EXP.Properties.VariableNames = lower(data_EXP.Properties.VariableNames);
data_EXP.rat = upper(strtrim(data_EXP.rat));  % uppercase for merges
if ~isempty(groupsToInclude_EXP)
    data_EXP = data_EXP(ismember(data_EXP.group, groupsToInclude_EXP), :);
end
nSessions_EXP = height(data_EXP);

% Initialize arrays for explore metrics
exploreLatency = NaN(nSessions_EXP,1);
explore_abs_decrease = NaN(nSessions_EXP,1);
explore_pct_decrease = NaN(nSessions_EXP,1);
explore_candidate_slope = NaN(nSessions_EXP,1);

for i = 1:nSessions_EXP
    reversalTrial = data_EXP.reversal_trial(i);
    
    % Get CPs from before and after reversal
    cp_before_trials = [data_EXP.cp1_trial_before(i), data_EXP.cp2_trial_before(i), ...
                        data_EXP.cp3_trial_before(i), data_EXP.cp4_trial_before(i)];
    cp_before_slopes = [data_EXP.cp1_slope_before(i), data_EXP.cp2_slope_before(i), ...
                        data_EXP.cp3_slope_before(i), data_EXP.cp4_slope_before(i)];
    cp_after_trials = [data_EXP.cp1_trial_after(i), data_EXP.cp2_trial_after(i), ...
                       data_EXP.cp3_trial_after(i), data_EXP.cp4_trial_after(i)];
    cp_after_slopes = [data_EXP.cp1_slope_after(i), data_EXP.cp2_slope_after(i), ...
                       data_EXP.cp3_slope_after(i), data_EXP.cp4_slope_after(i)];
    
    % Combine CPs into one list
    cp_all_trials = [cp_before_trials, cp_after_trials];
    cp_all_slopes = [cp_before_slopes, cp_after_slopes];
    
    % Remove NaN entries
    validIdx = ~isnan(cp_all_trials) & ~isnan(cp_all_slopes);
    cp_all_trials = cp_all_trials(validIdx);
    cp_all_slopes = cp_all_slopes(validIdx);
    
    % Sort CPs in ascending order by trial number
    [sortedTrials, sortIdx] = sort(cp_all_trials);
    sortedSlopes = cp_all_slopes(sortIdx);
    
    % Restrict to CPs occurring at or after (reversalTrial - maxTrialsBeforeReversal_EXP)
    window_start = reversalTrial - maxTrialsBeforeReversal_EXP;
    windowIdx = sortedTrials >= window_start;
    sortedTrials = sortedTrials(windowIdx);
    sortedSlopes = sortedSlopes(windowIdx);
    
    % Initialize candidate metrics as NaN
    candidateLatency = NaN;
    candidateAbsDecrease = NaN;
    candidate_pct_decrease = NaN;
    candidateSlope = NaN;
    
    % Look for the first decrease in slope (from the second CP onward)
    for j = 2:length(sortedTrials)
        if sortedSlopes(j) < sortedSlopes(j-1)
            % If a threshold is provided, enforce it
            if ~isempty(exploreSlopeThreshold_EXP) && ~isinf(exploreSlopeThreshold_EXP)
                if sortedSlopes(j) >= exploreSlopeThreshold_EXP
                    continue;
                end
            end
            % Candidate CP found; compute metrics
            candidateLatency = sortedTrials(j) - reversalTrial;  % (EPF1)
            candidateAbsDecrease = sortedSlopes(j-1) - sortedSlopes(j);
            if sortedSlopes(j-1) ~= 0
                candidate_pct_decrease = candidateAbsDecrease / sortedSlopes(j-1) * 100;
            end
            candidateSlope = sortedSlopes(j);
            break;
        end
    end
    
    exploreLatency(i) = candidateLatency;
    explore_abs_decrease(i) = candidateAbsDecrease;
    explore_pct_decrease(i) = candidate_pct_decrease;
    explore_candidate_slope(i) = candidateSlope;
end

% Append explore measures to the explore data table
data_EXP.explore_latency = exploreLatency;
data_EXP.explore_abs_decrease = explore_abs_decrease;
data_EXP.explore_pct_decrease = explore_pct_decrease;
data_EXP.explore_candidate_slope = explore_candidate_slope;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Aggregate to Rat-Level with Outlier Filtering and Pivot to Wide Format

% --- LF Measures Aggregation ---
[G, ratID, ratAge, ratGroup] = findgroups(data_LF.rat, data_LF.age, data_LF.group);
ratLF1 = splitapply(@(x) mean(x(~isoutlier(x)), 'omitnan'), data_LF.cp_before_latency, G);
T_LF1 = table(ratID, ratAge, ratLF1);
T_wide_LF1 = unstack(T_LF1, 'ratLF1', 'ratAge');
varNames = T_wide_LF1.Properties.VariableNames; varNames{1} = 'rat';
for k = 2:length(varNames)
    ageStr = varNames{k};
    if startsWith(ageStr, 'x'), ageStr = extractAfter(ageStr, 1); end
    varNames{k} = ['LF1_P' ageStr];
end
T_wide_LF1.Properties.VariableNames = varNames;

ratLF2_1 = splitapply(@(x) mean(x(~isoutlier(x)), 'omitnan'), data_LF.raw_first_candidate_slope, G);
T_LF2_1 = table(ratID, ratAge, ratLF2_1);
T_wide_LF2_1 = unstack(T_LF2_1, 'ratLF2_1', 'ratAge');
varNames = T_wide_LF2_1.Properties.VariableNames; varNames{1} = 'rat';
for k = 2:length(varNames)
    ageStr = varNames{k};
    if startsWith(ageStr, 'x'), ageStr = extractAfter(ageStr, 1); end
    varNames{k} = ['LF2_1_P' ageStr];
end
T_wide_LF2_1.Properties.VariableNames = varNames;

ratLF2_2 = splitapply(@(x) mean(x(~isoutlier(x)), 'omitnan'), data_LF.weighted_candidate_slope, G);
T_LF2_2 = table(ratID, ratAge, ratLF2_2);
T_wide_LF2_2 = unstack(T_LF2_2, 'ratLF2_2', 'ratAge');
varNames = T_wide_LF2_2.Properties.VariableNames; varNames{1} = 'rat';
for k = 2:length(varNames)
    ageStr = varNames{k};
    if startsWith(ageStr, 'x'), ageStr = extractAfter(ageStr, 1); end
    varNames{k} = ['LF2_2_P' ageStr];
end
T_wide_LF2_2.Properties.VariableNames = varNames;

ratLF2_3 = splitapply(@(x) mean(x(~isoutlier(x)), 'omitnan'), data_LF.weighted_candidate_slope_fig3, G);
T_LF2_3 = table(ratID, ratAge, ratLF2_3);
T_wide_LF2_3 = unstack(T_LF2_3, 'ratLF2_3', 'ratAge');
varNames = T_wide_LF2_3.Properties.VariableNames; varNames{1} = 'rat';
for k = 2:length(varNames)
    ageStr = varNames{k};
    if startsWith(ageStr, 'x'), ageStr = extractAfter(ageStr, 1); end
    varNames{k} = ['LF2_3_P' ageStr];
end
T_wide_LF2_3.Properties.VariableNames = varNames;

% --- ETF Measures Aggregation ---
[G, ratID, ratAge, ratGroup] = findgroups(data_ETF.rat, data_ETF.age, data_ETF.group);
ratETF1 = splitapply(@(x) mean(x(~isoutlier(x)), 'omitnan'), data_ETF.cp_before_latency_ETF, G);
T_ETF1 = table(ratID, ratAge, ratETF1);
T_wide_ETF1 = unstack(T_ETF1, 'ratETF1', 'ratAge');
varNames = T_wide_ETF1.Properties.VariableNames; varNames{1} = 'rat';
for k = 2:length(varNames)
    ageStr = varNames{k};
    if startsWith(ageStr, 'x'), ageStr = extractAfter(ageStr, 1); end
    varNames{k} = ['ETF1_P' ageStr];
end
T_wide_ETF1.Properties.VariableNames = varNames;

ratETF2_1 = splitapply(@(x) mean(x(~isoutlier(x)), 'omitnan'), data_ETF.raw_first_candidate_slope_ETF, G);
T_ETF2_1 = table(ratID, ratAge, ratETF2_1);
T_wide_ETF2_1 = unstack(T_ETF2_1, 'ratETF2_1', 'ratAge');
varNames = T_wide_ETF2_1.Properties.VariableNames; varNames{1} = 'rat';
for k = 2:length(varNames)
    ageStr = varNames{k};
    if startsWith(ageStr, 'x'), ageStr = extractAfter(ageStr, 1); end
    varNames{k} = ['ETF2_1_P' ageStr];
end
T_wide_ETF2_1.Properties.VariableNames = varNames;

ratETF2_2 = splitapply(@(x) mean(x(~isoutlier(x)), 'omitnan'), data_ETF.weighted_candidate_slope_ETF, G);
T_ETF2_2 = table(ratID, ratAge, ratETF2_2);
T_wide_ETF2_2 = unstack(T_ETF2_2, 'ratETF2_2', 'ratAge');
varNames = T_wide_ETF2_2.Properties.VariableNames; varNames{1} = 'rat';
for k = 2:length(varNames)
    ageStr = varNames{k};
    if startsWith(ageStr, 'x'), ageStr = extractAfter(ageStr, 1); end
    varNames{k} = ['ETF2_2_P' ageStr];
end
T_wide_ETF2_2.Properties.VariableNames = varNames;

ratETF2_3 = splitapply(@(x) mean(x(~isoutlier(x)), 'omitnan'), data_ETF.weighted_candidate_slope_fig3_ETF, G);
T_ETF2_3 = table(ratID, ratAge, ratETF2_3);
T_wide_ETF2_3 = unstack(T_ETF2_3, 'ratETF2_3', 'ratAge');
varNames = T_wide_ETF2_3.Properties.VariableNames; varNames{1} = 'rat';
for k = 2:length(varNames)
    ageStr = varNames{k};
    if startsWith(ageStr, 'x'), ageStr = extractAfter(ageStr, 1); end
    varNames{k} = ['ETF2_3_P' ageStr];
end
T_wide_ETF2_3.Properties.VariableNames = varNames;

% --- Explore Measures Aggregation ---
[G, ratID, ratAge, ratGroup] = findgroups(data_EXP.rat, data_EXP.age, data_EXP.group);
ratEPF1 = splitapply(@(x) mean(x(~isoutlier(x)), 'omitnan'), data_EXP.explore_latency, G);
T_EP_lat = table(ratID, ratAge, ratGroup, ratEPF1);
T_wide_EP_lat = unstack(T_EP_lat, 'ratEPF1', 'ratAge');
varNames = T_wide_EP_lat.Properties.VariableNames; varNames{1} = 'rat';
for k = 2:length(varNames)
    ageStr = varNames{k};
    if startsWith(ageStr, 'x'), ageStr = extractAfter(ageStr, 1); end
    varNames{k} = ['EPF1_P' ageStr];
end
T_wide_EP_lat.Properties.VariableNames = varNames;

ratEPF2_1 = splitapply(@(x) mean(x(~isoutlier(x)), 'omitnan'), data_EXP.explore_abs_decrease, G);
T_EP_2_1 = table(ratID, ratAge, ratEPF2_1);
T_wide_EP_2_1 = unstack(T_EP_2_1, 'ratEPF2_1', 'ratAge');
varNames = T_wide_EP_2_1.Properties.VariableNames; varNames{1} = 'rat';
for k = 2:length(varNames)
    ageStr = varNames{k};
    if startsWith(ageStr, 'x'), ageStr = extractAfter(ageStr, 1); end
    varNames{k} = ['EPF2_1_P' ageStr];
end
T_wide_EP_2_1.Properties.VariableNames = varNames;

ratEPF2_2 = splitapply(@(x) mean(x(~isoutlier(x)), 'omitnan'), data_EXP.explore_pct_decrease, G);
T_EP_2_2 = table(ratID, ratAge, ratEPF2_2);
T_wide_EP_2_2 = unstack(T_EP_2_2, 'ratEPF2_2', 'ratAge');
varNames = T_wide_EP_2_2.Properties.VariableNames; varNames{1} = 'rat';
for k = 2:length(varNames)
    ageStr = varNames{k};
    if startsWith(ageStr, 'x'), ageStr = extractAfter(ageStr, 1); end
    varNames{k} = ['EPF2_2_P' ageStr];
end
T_wide_EP_2_2.Properties.VariableNames = varNames;

ratEPF2_3 = splitapply(@(x) mean(x(~isoutlier(x)), 'omitnan'), data_EXP.explore_candidate_slope, G);
T_EP_2_3 = table(ratID, ratAge, ratEPF2_3);
T_wide_EP_2_3 = unstack(T_EP_2_3, 'ratEPF2_3', 'ratAge');
varNames = T_wide_EP_2_3.Properties.VariableNames; varNames{1} = 'rat';
for k = 2:length(varNames)
    ageStr = varNames{k};
    if startsWith(ageStr, 'x'), ageStr = extractAfter(ageStr, 1); end
    varNames{k} = ['EPF2_3_P' ageStr];
end
T_wide_EP_2_3.Properties.VariableNames = varNames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Rat Metadata
if ismember('sex', data_all.Properties.VariableNames) && ismember('group', data_all.Properties.VariableNames)
    T_meta = unique(data_all(:, {'rat', 'sex', 'group'}));
else
    T_meta = table(unique(data_all.rat), 'VariableNames', {'rat'});
end

% Already uppercase from data_all; ensure trimming
T_meta.rat = upper(strtrim(T_meta.rat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Merge All Wide Tables with Metadata and Add Reversal Measures

% Define merge options (without ConflictVariableNames)
opts = {'Keys','rat','MergeKeys',true,'Type','left'};

T_temp = outerjoin(T_meta, T_wide_LF1, opts{:});
T_temp = outerjoin(T_temp, T_wide_LF2_1, opts{:});
T_temp = outerjoin(T_temp, T_wide_LF2_2, opts{:});
T_temp = outerjoin(T_temp, T_wide_LF2_3, opts{:});
T_temp = outerjoin(T_temp, T_wide_ETF1, opts{:});
T_temp = outerjoin(T_temp, T_wide_ETF2_1, opts{:});
T_temp = outerjoin(T_temp, T_wide_ETF2_2, opts{:});
T_temp = outerjoin(T_temp, T_wide_ETF2_3, opts{:});
T_temp = outerjoin(T_temp, T_wide_EP_lat, opts{:});
T_temp = outerjoin(T_temp, T_wide_EP_2_1, opts{:});
T_temp = outerjoin(T_temp, T_wide_EP_2_2, opts{:});
T_temp = outerjoin(T_temp, T_wide_EP_2_3, opts{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load and Prepare Reversal Data
revData = readtable('/Users/olly/Desktop/ASI/ASI/MAT Files/avg_reversals_per_rat.csv');
revData.Properties.VariableNames = lower(revData.Properties.VariableNames);

% Convert rat IDs to uppercase for merging
revData.rat = upper(strtrim(revData.rat));

% --- Filter Reversal Data ---
excluded_rats = {'DEV019020016','DEV005017010'};  % uppercase to match
included_groups = groupsToInclude;  % e.g., {'SOC'}
included_sexes = {'M','F'};           % options: {'M','F'}
included_ages = [30 50 70 90 120 190];

% Remove excluded rats
revData(ismember(revData.rat, excluded_rats), :) = [];
% Filter by group
revData = revData(ismember(revData.group, included_groups), :);

% Filter by sex using core data
core = readtable('/Users/olly/Desktop/ASI/ASI/MAT Files/core_ASI_data.csv');
core.Properties.VariableNames = lower(core.Properties.VariableNames);
core.rat = upper(strtrim(core.rat));
sex_lookup = containers.Map(core.rat, core.sex);
keep_mask = false(height(revData), 1);
for i = 1:height(revData)
    rat_id = revData.rat{i};
    if isKey(sex_lookup, rat_id)
        keep_mask(i) = ismember(sex_lookup(rat_id), included_sexes);
    end
end
revData = revData(keep_mask, :);

% Convert reversal measure columns (except 'rat' and 'group') to numeric if needed
revVarNames = revData.Properties.VariableNames;
for k = 1:length(revVarNames)
    if ~ismember(revVarNames{k}, {'rat', 'group'})
        if ~isnumeric(revData.(revVarNames{k}))
            revData.(revVarNames{k}) = str2double(revData.(revVarNames{k}));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Final Merge with Reversal Data
T_final = outerjoin(T_temp, revData, 'Keys','rat','MergeKeys',true,'Type','left');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Export Final Table to CSV
exportFileName = sprintf('data_test_%.1f.csv', critVal);
exportFilePath = fullfile('/Users/olly/Desktop/ASI/ASI/MAT Files/', exportFileName);
writetable(T_final, exportFilePath);

fprintf('\nExported combined regression data to:\n%s\n\n', exportFilePath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local helper function for computing mean ignoring NaNs
function y = localNanMean(x, dim)
    if nargin < 2
        dim = find(size(x)~=1, 1);
        if isempty(dim), dim = 1; end
    end
    sumX = sum(x, dim, 'omitnan');
    countX = sum(~isnan(x), dim);
    y = sumX ./ countX;
    y(countX == 0) = NaN;
end