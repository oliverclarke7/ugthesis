%% investigatingsinglerat_with_adjustables.m
% This script loads your data, applies global adjustable filters, computes
% the CP_before latency using specified thresholds, and then queries the data
% for a particular rat and age. Optionally, you can further filter by reversal
% name (inclusion and exclusion) at the query level.
%
% Adjustable parameters for CP_before latency and reversal filtering are defined
% below. You can modify these values as needed.

clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adjustable Parameters

% Parameters for CP_before latency computation:
critVal = 1.5;                % Adjustable critical value
slopeThreshold = 0.6;         % Slope threshold for candidate CP
minTrialsBeforeNextCP = 0;    % Minimum trials gap to the next CP (adjustable)
minTrialsBeforeReversal = 0;  % Minimum trials gap to the reversal (adjustable)
maxTrialsBeforeReversal = [35]; % Maximum trials before reversal (set [] or Inf for no limit)

% Global reversal filter parameters:
includedReversal = {};         % Specify reversals to include (e.g., {'REV3_REV4'}); leave empty to ignore
excludedReversal = 'REV1_REV2';  % Exclude this reversal_name if inclusion filter is empty

% Group filter:
groupsToInclude = {'SOC'};     % Only include sessions from group 'SOC' (modify as needed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Query Parameters (for a specific rat/age and further reversal filtering)
targetRat = 'DEV005017009';            % Specify the rat (as a string)
targetAge = 120;               % Specify the age (e.g., 120)
targetReversalQuery = '';      % Specify a reversal name to include in the query (e.g., 'REV3_REV4'); leave empty for no additional filter
targetReversalExcludeQuery = {}; % Specify reversal types to exclude in the query (e.g., {'REV1_REV2'}); leave empty if none

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data File Settings
dataDir = '/Users/olly/Desktop/ASI/ASI/MAT Files/';
dataFile = sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal);
filePath = fullfile(dataDir, dataFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data and Standardize Variable Names
data = readtable(filePath);
data.Properties.VariableNames = lower(data.Properties.VariableNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply Global Reversal and Group Filters

% Global reversal filter: if an inclusion filter is provided, use it; otherwise, apply exclusion.
if ~isempty(includedReversal)
    data = data(ismember(data.reversal_name, includedReversal), :);
elseif ~isempty(excludedReversal)
    data = data(~ismember(data.reversal_name, excludedReversal), :);
end

% Filter by group.
data = data(ismember(data.group, groupsToInclude), :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute CP_before Latency for Each Session Using Adjustable Parameters

cpBeforeLatency = NaN(height(data), 1);

for i = 1:height(data)
    % Extract slopes and trial numbers for CPs before reversal
    slopes_before = [data.cp1_slope_before(i), data.cp2_slope_before(i), ...
                     data.cp3_slope_before(i), data.cp4_slope_before(i)];
    trials_before = [data.cp1_trial_before(i), data.cp2_trial_before(i), ...
                     data.cp3_trial_before(i), data.cp4_trial_before(i)];
    reversalTrial = data.reversal_trial(i);
    
    % Identify candidate indices where slope >= slopeThreshold (ignoring NaNs)
    candidateIndices = find(slopes_before >= slopeThreshold & ~isnan(slopes_before));
    
    for j = 1:length(candidateIndices)
        idx = candidateIndices(j);
        candidateTrial = trials_before(idx);
        gapToReversal = reversalTrial - candidateTrial;
        
        % Condition 1: Must occur at least minTrialsBeforeReversal trials before reversal.
        if gapToReversal < minTrialsBeforeReversal
            continue;
        end
        
        % Condition 2: If maxTrialsBeforeReversal is set, candidate must occur no more than that many trials before reversal.
        if ~isempty(maxTrialsBeforeReversal) && ~isinf(maxTrialsBeforeReversal)
            if gapToReversal > maxTrialsBeforeReversal
                continue;
            end
        end
        
        % Condition 3: Check gap to next CP (if available) or to reversal.
        if idx < length(trials_before) && ~isnan(trials_before(idx+1))
            nextEventTrial = trials_before(idx+1);
        else
            nextEventTrial = reversalTrial;
        end
        
        if (nextEventTrial - candidateTrial) < minTrialsBeforeNextCP
            continue;
        end
        
        % Candidate meets all conditions: record its trial number.
        cpBeforeLatency(i) = candidateTrial;
        break; % Only take the first candidate
    end
end

% Add computed latency to the data table.
data.cp_before_latency = cpBeforeLatency;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Query the Data for a Specific Rat, Age, and Optional Reversal Filtering

% First, filter by rat and age.
idxQuery = strcmp(data.rat, targetRat) & (data.age == targetAge);

% Then apply additional query filters if provided.
if ~isempty(targetReversalQuery)
    idxQuery = idxQuery & strcmp(data.reversal_name, targetReversalQuery);
end
if ~isempty(targetReversalExcludeQuery)
    idxQuery = idxQuery & ~ismember(data.reversal_name, targetReversalExcludeQuery);
end

queryData = data(idxQuery, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display the Results

% Determine which columns to display (including cp_before_latency if available)
if ismember('cp_before_latency', data.Properties.VariableNames)
    displayCols = {'rat', 'age', 'reversal_name', 'cp_before_latency', 'reversal_trial'};
else
    displayCols = {'rat', 'age', 'reversal_name', 'reversal_trial'};
end

if isempty(queryData)
    fprintf('No data found for rat %s, age %d.\n', targetRat, targetAge);
else
    fprintf('Data for rat %s, age %d:\n', targetRat, targetAge);
    disp(queryData(:, displayCols));
end