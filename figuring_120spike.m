%% exploit_find_1.m
% This script computes CP_before latency for candidate CPs in EXPLOIT sessions.
% For EXPLOIT, you can now either specify a reversal inclusion filter or an exclusion filter.
% After computing session‐level CP_before latencies, we average within each rat (by age)
% and then average across rats (by age and group).
%
% Additionally, we compute the reversal trial number using the same within‐rat averaging 
% procedure, and then aggregate across rats. Both metrics are plotted with SEM error bars.

%% PARAMETERS
critVal = 1.5;                % Adjustable critical value
slopeThreshold = 0.6;         % Slope threshold for candidate CP
minTrialsBeforeNextCP = 0;    % Minimum trials gap to the next CP (adjustable)
minTrialsBeforeReversal = 0;  % Minimum trials gap to the reversal (adjustable)
maxTrialsBeforeReversal = [35]; % Maximum trials before reversal (set [] or Inf for no limit)

% Reversal filter parameters:
includedReversal = {};         % Specify reversals to include (e.g., {'REV3_REV4'}); leave empty to ignore
excludedReversal = 'REV1_REV2';  % Exclude this reversal_name if inclusion filter is empty

groupsToInclude = {'SOC'};     % Only include sessions from group 'SOC' (modify as needed)

dataDir = '/Users/olly/Desktop/ASI/ASI/MAT Files/';
dataFile = sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal);
filePath = fullfile(dataDir, dataFile);

%% Load data
data = readtable(filePath);
data.Properties.VariableNames = lower(data.Properties.VariableNames);

%% Filter data for EXPLOIT
% Apply reversal inclusion/exclusion logic:
if ~isempty(includedReversal)
    data = data(ismember(data.reversal_name, includedReversal), :);
elseif ~isempty(excludedReversal)
    data = data(~ismember(data.reversal_name, excludedReversal), :);
end

% Keep only rows whose group is in groupsToInclude.
data = data(ismember(data.group, groupsToInclude), :);

% Extract unique ages and groups (for later plotting)
ages = unique(data.age);
groups = unique(data.group);

%% Compute CP_before latency for candidate CPs meeting all criteria (session-level)
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

% Add computed latency to the table.
data.cp_before_latency = cpBeforeLatency;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Within-Rat Averaging: Compute rat-level averages for CP_before latency
% (Assumes that the data table includes a "rat" column for subject identification.)
[G, ratID, ratAge, ratGroup] = findgroups(data.rat, data.age, data.group);
ratLatency = splitapply(@(x) mean(x, 'omitnan'), data.cp_before_latency, G);

% Also compute within-rat average reversal trial number.
ratReversal = splitapply(@(x) mean(x, 'omitnan'), data.reversal_trial, G);

% Create a table of rat-level averages.
T_rat = table(ratID, ratAge, ratGroup, ratLatency, ratReversal);

%% Display rat-level averages for rats at age 120 (P120)
idx_P120 = (T_rat.ratAge == 120);
ratP120 = T_rat(idx_P120, :);
disp('Rat-level average latency at P120:');
disp(ratP120);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Across-Rat Aggregation by Age and Group for CP_before latency
uniqueAges = unique(T_rat.ratAge);
uniqueGroups = unique(T_rat.ratGroup);
nAges = numel(uniqueAges);
nGroups = numel(uniqueGroups);

meanLatency = NaN(nAges, nGroups);
stderrLatency = NaN(nAges, nGroups);

for g = 1:nGroups
    for a = 1:nAges
        idx = (T_rat.ratAge == uniqueAges(a)) & strcmp(T_rat.ratGroup, uniqueGroups{g});
        if sum(idx) > 0
            meanLatency(a, g) = mean(T_rat.ratLatency(idx), 'omitnan');
            stderrLatency(a, g) = std(T_rat.ratLatency(idx), 'omitnan') / sqrt(sum(idx));
        else
            meanLatency(a, g) = NaN;
            stderrLatency(a, g) = NaN;
        end
    end
end

% Overall CP_before latency (averaged across groups)
overallMeanLatency = localNanMean(meanLatency, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Across-Rat Aggregation by Age and Group for reversal trial numbers
meanReversalRat = NaN(nAges, nGroups);
stderrReversalRat = NaN(nAges, nGroups);

for g = 1:nGroups
    for a = 1:nAges
        idx = (T_rat.ratAge == uniqueAges(a)) & strcmp(T_rat.ratGroup, uniqueGroups{g});
        if sum(idx) > 0
            meanReversalRat(a, g) = mean(T_rat.ratReversal(idx), 'omitnan');
            stderrReversalRat(a, g) = std(T_rat.ratReversal(idx), 'omitnan') / sqrt(sum(idx));
        else
            meanReversalRat(a, g) = NaN;
            stderrReversalRat(a, g) = NaN;
        end
    end
end

% Overall reversal trial number (averaged across groups)
overallMeanReversal = localNanMean(meanReversalRat, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot: CP Before Latency by Age and Group with SEM Error Bars,
%%       and reversal trial number (averaged in the same way)
figure;
hold on;
% Plot each group's CP_before latency with error bars.
for g = 1:nGroups
    errorbar(uniqueAges, meanLatency(:, g), stderrLatency(:, g), '-o', 'LineWidth', 1.5);
end
% Plot overall CP_before latency (black solid line).
errorbar(uniqueAges, overallMeanLatency, zeros(size(overallMeanLatency)), 'k-', 'LineWidth', 2);

% Now plot reversal trial number averaged by rat as a red dashed line with error bars.
for g = 1:nGroups
    errorbar(uniqueAges, meanReversalRat(:, g), stderrReversalRat(:, g), 'r--', 'LineWidth', 1.5);
end
errorbar(uniqueAges, overallMeanReversal, zeros(size(overallMeanReversal)), 'r-', 'LineWidth', 2);

% Build dynamic legend entries and title.
if ~isempty(includedReversal)
    if iscell(includedReversal)
        incStr = strjoin(includedReversal, ', ');
    else
        incStr = includedReversal;
    end
    revStr = ['Including ', incStr];
elseif ~isempty(excludedReversal)
    if iscell(excludedReversal)
        excStr = strjoin(excludedReversal, ', ');
    else
        excStr = excludedReversal;
    end
    revStr = ['Excluding ', excStr];
else
    revStr = 'All Reversals';
end

legendEntries = [uniqueGroups; {'Overall'}; uniqueGroups; {'Overall Reversal'}];
legend(legendEntries, 'Location', 'best');
title(sprintf('CP Before Latency (Slope ≥ %.1f, Crit %.1f)\n(EXPLOIT: %s, Groups=%s)', ...
    slopeThreshold, critVal, revStr, strjoin(groupsToInclude, ',')));
xlabel('Age');
ylabel('Trial Number (Latency relative to Trial 0) / Reversal Trial Number');
hold off;

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