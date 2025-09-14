%% export_ETF2_metrics_singleRowPerRat.m
% This script:
% 1) Loads and filters session data for EXPLOIT.
% 2) Computes three CP metrics for each session:
%    - Raw first candidate slope (Graph 1)
%    - Weighted candidate slope (Graph 2)
%    - Weighted slope for all CPs from first candidate to reversal (Graph 3)
% 3) Averages those metrics by (rat, age).
% 4) Pivots so that each rat appears exactly once (one row per rat),
%    with columns for each metric at each age (e.g., ETF2_1_P30, ETF2_2_P30, ETF2_3_P30).
% 5) Writes the result to 'ETF2_test.csv' in the specified folder.
%
% Adjustable parameters are at the top of the script.

clear; clc;

%% PARAMETERS
critVal = 3.2;                % Critical value for file selection
slopeThreshold = 0.6;         % Minimum slope for candidate CP
minTrialsBeforeNextCP = 5;    % Minimum trial gap to next CP
minTrialsBeforeReversal = 0;  % Minimum trials gap to the reversal
maxTrialsBeforeReversal = 35; % Maximum trials before reversal (set [] or Inf for no limit)

% Filtering parameters:
excludedReversal = 'REV1_REV2';  % Exclude sessions with this reversal name
groupsToInclude = {'SOC'};       % Include only these groups

% File paths
dataDir = '/Users/olly/Desktop/ASI/ASI/MAT Files/';
dataFile = sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal);
filePath = fullfile(dataDir, dataFile);

% Output file
outputFile = fullfile(dataDir, 'ETF2_test.csv');

%% 1) LOAD DATA
data = readtable(filePath);
data.Properties.VariableNames = lower(data.Properties.VariableNames);

%% 2) FILTER DATA (EXPLOIT ONLY, EXCLUDING SPECIFIC REVERSAL AND GROUPS)
data = data(~strcmp(data.reversal_name, excludedReversal), :);      % Exclude specified reversal
data = data(ismember(data.group, groupsToInclude), :);             % Keep only specified group(s)

%% 3) COMPUTE SESSION-LEVEL METRICS

nSessions = height(data);

% Preallocate
rawFirstCandidateSlope = NaN(nSessions, 1);        % Graph 1
weightedCandidateSlope = NaN(nSessions, 1);        % Graph 2
weightedCandidateSlopeFig3 = NaN(nSessions, 1);    % Graph 3

% --- Loop for Graphs 1 & 2 ---
for i = 1:nSessions
    slopes_before = [data.cp1_slope_before(i), data.cp2_slope_before(i), ...
                     data.cp3_slope_before(i), data.cp4_slope_before(i)];
    trials_before = [data.cp1_trial_before(i), data.cp2_trial_before(i), ...
                     data.cp3_trial_before(i), data.cp4_trial_before(i)];
    reversalTrial = data.reversal_trial(i);
    
    candidateSlopes = [];
    candidateWeights = [];
    firstFound = false;
    
    for j = 1:length(slopes_before)
        currSlope = slopes_before(j);
        currTrial = trials_before(j);
        if isnan(currSlope) || isnan(currTrial)
            continue;
        end
        if currSlope < slopeThreshold
            continue;
        end
        gapToReversal = reversalTrial - currTrial;
        if gapToReversal < minTrialsBeforeReversal
            continue;
        end
        if ~isempty(maxTrialsBeforeReversal) && ~isinf(maxTrialsBeforeReversal)
            if gapToReversal > maxTrialsBeforeReversal
                continue;
            end
        end
        % Next event is the next CP or the reversal
        if j < length(trials_before) && ~isnan(trials_before(j+1))
            nextEventTrial = trials_before(j+1);
        else
            nextEventTrial = reversalTrial;
        end
        if (nextEventTrial - currTrial) < minTrialsBeforeNextCP
            continue;
        end
        
        % It's a valid candidate
        candidateSlopes(end+1) = currSlope; %#ok<AGROW>
        candidateWeights(end+1) = nextEventTrial - currTrial; %#ok<AGROW>
        
        if ~firstFound
            rawFirstCandidateSlope(i) = currSlope;  % Graph 1
            firstFound = true;
        end
    end
    
    % Graph 2: Weighted average if multiple candidates
    if ~isempty(candidateSlopes)
        weightedCandidateSlope(i) = sum(candidateSlopes .* candidateWeights) / sum(candidateWeights);
    end
end

% --- Loop for Graph 3 (weighted slope for all CPs from first candidate to reversal) ---
for i = 1:nSessions
    slopes_before = [data.cp1_slope_before(i), data.cp2_slope_before(i), ...
                     data.cp3_slope_before(i), data.cp4_slope_before(i)];
    trials_before = [data.cp1_trial_before(i), data.cp2_trial_before(i), ...
                     data.cp3_trial_before(i), data.cp4_trial_before(i)];
    reversalTrial = data.reversal_trial(i);
    
    firstCandidateIndex = [];
    for j = 1:length(slopes_before)
        currSlope = slopes_before(j);
        currTrial = trials_before(j);
        if isnan(currSlope) || isnan(currTrial)
            continue;
        end
        if currSlope < slopeThreshold
            continue;
        end
        gapToReversal = reversalTrial - currTrial;
        if gapToReversal < minTrialsBeforeReversal
            continue;
        end
        if ~isempty(maxTrialsBeforeReversal) && ~isinf(maxTrialsBeforeReversal)
            if gapToReversal > maxTrialsBeforeReversal
                continue;
            end
        end
        if j < length(trials_before) && ~isnan(trials_before(j+1))
            nextEventTrial = trials_before(j+1);
        else
            nextEventTrial = reversalTrial;
        end
        if (nextEventTrial - currTrial) < minTrialsBeforeNextCP
            continue;
        end
        firstCandidateIndex = j;
        break;
    end
    
    if ~isempty(firstCandidateIndex)
        candidateSlopesFig3 = [];
        candidateWeightsFig3 = [];
        for j = firstCandidateIndex:length(slopes_before)
            if isnan(trials_before(j))
                continue;
            end
            if trials_before(j) >= reversalTrial
                continue;
            end
            slope_j = slopes_before(j);
            if (j < length(trials_before)) && (~isnan(trials_before(j+1))) && (trials_before(j+1) < reversalTrial)
                weight_j = trials_before(j+1) - trials_before(j);
            else
                weight_j = reversalTrial - trials_before(j);
            end
            candidateSlopesFig3(end+1) = slope_j;      %#ok<AGROW>
            candidateWeightsFig3(end+1) = weight_j;    %#ok<AGROW>
        end
        if ~isempty(candidateSlopesFig3)
            weightedCandidateSlopeFig3(i) = sum(candidateSlopesFig3 .* candidateWeightsFig3) / sum(candidateWeightsFig3);
        end
    end
end

% Store results back into data
data.raw_first_candidate_slope = rawFirstCandidateSlope;
data.weighted_candidate_slope = weightedCandidateSlope;
data.weighted_candidate_slope_fig3 = weightedCandidateSlopeFig3;

%% 4) AVERAGE ACROSS SESSIONS AT (RAT, AGE) LEVEL

% Group by rat and age (ignoring group to ensure one row per rat)
[G, ratID, ratAge] = findgroups(data.rat, data.age);

% Compute means for each group
ratRawSlope = splitapply(@(x) mean(x, 'omitnan'), data.raw_first_candidate_slope, G);
ratWeightedSlope = splitapply(@(x) mean(x, 'omitnan'), data.weighted_candidate_slope, G);
ratWeightedSlopeFig3 = splitapply(@(x) mean(x, 'omitnan'), data.weighted_candidate_slope_fig3, G);

% Build a table of unique (rat, age) combos
T_ratAge = table(ratID, ratAge, ratRawSlope, ratWeightedSlope, ratWeightedSlopeFig3);

%% 5) PIVOT SO THAT EACH RAT IS A SINGLE ROW
uniqueRats = unique(T_ratAge.ratID, 'stable');  % keep stable ordering
uniqueAges = unique(T_ratAge.ratAge);
nRats = numel(uniqueRats);
nAges = numel(uniqueAges);

% Initialize final arrays
ETF2_1 = NaN(nRats, nAges);  % raw slope
ETF2_2 = NaN(nRats, nAges);  % weighted slope
ETF2_3 = NaN(nRats, nAges);  % weighted slope (all CPs)

% Fill these arrays for each rat.
for iRat = 1:nRats
    currRat = uniqueRats{iRat};  % Note: currRat is a string
    % Use strcmp to compare strings
    idxRat = strcmp(T_ratAge.ratID, currRat);
    % For each age the rat has data for
    rowsForRat = T_ratAge(idxRat, :);
    for k = 1:height(rowsForRat)
        thisAge = rowsForRat.ratAge(k);
        ageCol = find(uniqueAges == thisAge);
        
        ETF2_1(iRat, ageCol) = rowsForRat.ratRawSlope(k);
        ETF2_2(iRat, ageCol) = rowsForRat.ratWeightedSlope(k);
        ETF2_3(iRat, ageCol) = rowsForRat.ratWeightedSlopeFig3(k);
    end
end

% Create variable names for columns
colNames_1 = arrayfun(@(a) sprintf('ETF2_1_P%d', a), uniqueAges, 'UniformOutput', false);
colNames_2 = arrayfun(@(a) sprintf('ETF2_2_P%d', a), uniqueAges, 'UniformOutput', false);
colNames_3 = arrayfun(@(a) sprintf('ETF2_3_P%d', a), uniqueAges, 'UniformOutput', false);

% Convert to tables
T_ETF2_1 = array2table(ETF2_1, 'VariableNames', colNames_1);
T_ETF2_2 = array2table(ETF2_2, 'VariableNames', colNames_2);
T_ETF2_3 = array2table(ETF2_3, 'VariableNames', colNames_3);

% Final table: one row per rat, with columns for each metric/age combination
T_out = table(uniqueRats, 'VariableNames', {'rat'});
T_out = [T_out, T_ETF2_1, T_ETF2_2, T_ETF2_3];

%% 6) WRITE TO CSV
writetable(T_out, outputFile);
fprintf('Exported %d rows (one per rat) to %s\n', height(T_out), outputFile);