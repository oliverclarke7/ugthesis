%% Script: Add Rat-Level Average CP_before Latency Columns Using the Correct Data Sheet
% This script reads session-level data from the CSV file used in your original code,
% computes CP_before latency for candidate change points (CPs) in EXPLOIT sessions,
% averages the latency within each rat (grouped by age), and pivots the data to a
% wide-format table with one row per rat and columns for each age (e.g., ETF1_P30, ETF1_P50, etc.).
% It then saves the result to "ETF1_test.csv".

%% PARAMETERS
% Critical value and file I/O parameters (using the same sheet as in your original code)
critVal = 3.2;                % Adjustable critical value
dataDir = '/Users/olly/Desktop/ASI/ASI/MAT Files/';
dataFile = sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal);
inputFile = fullfile(dataDir, dataFile);
outputFile = fullfile(dataDir, 'ETF1_test.csv');  % Save as ETF1_test.csv

% CP_before latency criteria
slopeThreshold        = 0.6;   % Slope threshold for candidate CP
minTrialsBeforeNextCP = 0;     % Minimum trials gap to the next CP
minTrialsBeforeReversal = 0;   % Minimum trials gap to reversal
maxTrialsBeforeReversal = 35;  % Maximum allowed trials before reversal (set [] or Inf for no limit)

% Filtering parameters for EXPLOIT sessions:
excludedReversal = 'REV1_REV2'; % Exclude sessions with this reversal name
groupsToInclude  = {'SOC'};     % Only include sessions from group 'SOC'

%% Load Data
data = readtable(inputFile);
data.Properties.VariableNames = lower(data.Properties.VariableNames);

% Check for the reversal_name column (try alternative 'reversalname' if needed)
if ~ismember('reversal_name', data.Properties.VariableNames)
    if ismember('reversalname', data.Properties.VariableNames)
        data.reversal_name = data.reversalname;
    else
        error(['The CSV file does not contain a column named "reversal_name" or "reversalname". ' ...
               'Available columns: %s'], strjoin(data.Properties.VariableNames, ', '));
    end
end

%% Filter Data for EXPLOIT Sessions
% Exclude sessions with the specified reversal and include only designated groups.
data = data(~strcmp(data.reversal_name, excludedReversal), :);
data = data(ismember(data.group, groupsToInclude), :);

%% Compute CP_before Latency for Each Session
% Preallocate vector for CP_before latency.
cpBeforeLatency = NaN(height(data), 1);

for i = 1:height(data)
    % Extract slopes and trial numbers for candidate CPs before reversal.
    slopes_before = [data.cp1_slope_before(i), data.cp2_slope_before(i), ...
                     data.cp3_slope_before(i), data.cp4_slope_before(i)];
    trials_before = [data.cp1_trial_before(i), data.cp2_trial_before(i), ...
                     data.cp3_trial_before(i), data.cp4_trial_before(i)];
    reversalTrial = data.reversal_trial(i);
    
    % Identify candidate CP indices (where slope >= threshold and not NaN)
    candidateIndices = find(slopes_before >= slopeThreshold & ~isnan(slopes_before));
    
    for j = 1:length(candidateIndices)
        idx = candidateIndices(j);
        candidateTrial = trials_before(idx);
        gapToReversal = reversalTrial - candidateTrial;
        
        % Condition 1: Must occur at least minTrialsBeforeReversal trials before reversal.
        if gapToReversal < minTrialsBeforeReversal
            continue;
        end
        
        % Condition 2: Must occur no more than maxTrialsBeforeReversal trials before reversal.
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
        break;
    end
end
data.cp_before_latency = cpBeforeLatency;

%% Within-Rat Averaging: Compute Averages by Rat and Age
% Group data by rat, age, and group.
[G, ratID, ratAge, ratGroup] = findgroups(data.rat, data.age, data.group);
ratLatency = splitapply(@(x) mean(x, 'omitnan'), data.cp_before_latency, G);

% Create a table with rat-level averages (each row is one rat at a given age).
T_rat = table(ratID, ratAge, ratGroup, ratLatency);

%% Pivot to Wide-Format: One Row per Rat, Columns for Each Age
% unstack() will produce columns named x30, x50, etc. when 'ratAge' is numeric.
T_wide = unstack(T_rat, 'ratLatency', 'ratAge');

%% Rename columns x30, x50, etc. to ETF1_P30, ETF1_P50, etc.
varNames = T_wide.Properties.VariableNames;
for k = 1:numel(varNames)
    % If a column name starts with 'x' (e.g., x30), rename it.
    if startsWith(varNames{k}, 'x')
        ageStr = varNames{k}(2:end);  % e.g., '30' from 'x30'
        newName = ['ETF1_P' ageStr];  % e.g., 'ETF1_P30'
        T_wide = renamevars(T_wide, varNames{k}, newName);
    end
end

%% Save the Wide-Format Table to CSV
writetable(T_wide, outputFile);
disp('Wide-format rat-level data saved successfully as ETF1_test.csv.');