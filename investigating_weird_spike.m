%% investigating_weird_spike.m
% This script investigates an unexpected spike in reversal trial numbers
% for EXPLOIT sessions at a target age (e.g., Age = 120), using the same 
% within‐rat averaging procedure as in exploit_find_1.m.

%% PARAMETERS
critVal = 3.2;                      % Critical value (should match your file)
excludedReversal = 'REV1_REV2';     % Exclude sessions with this reversal_name (EXPLOIT filter)
groupsToInclude = {'SOC'};          % Groups of interest
suspiciousAge = 120;                % Target age for investigation

dataDir = '/Users/olly/Desktop/ASI/ASI/MAT Files/';
dataFile = sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal);
filePath = fullfile(dataDir, dataFile);

%% Load data
data = readtable(filePath);
data.Properties.VariableNames = lower(data.Properties.VariableNames);

%% Filter data for EXPLOIT sessions
% Exclude sessions with reversal_name equal to excludedReversal.
data = data(~strcmp(data.reversal_name, excludedReversal), :);
% Keep only rows whose group is in groupsToInclude.
data = data(ismember(data.group, groupsToInclude), :);

%% Identify sessions at the target age (Age = 120)
idxSuspicious = (data.age == suspiciousAge);

if ~any(idxSuspicious)
    fprintf('No sessions found at Age = %d.\n', suspiciousAge);
    return;
end

suspiciousData = data(idxSuspicious, :);

fprintf('\n=== Investigation of Age = %d ===\n', suspiciousAge);
fprintf('Number of sessions at this age: %d\n', height(suspiciousData));

% Identify unique rats present at this age.
uniqueRatsAtSuspiciousAge = unique(suspiciousData.rat);
fprintf('Number of unique rats at this age: %d\n', numel(uniqueRatsAtSuspiciousAge));
disp('Rats present:');
disp(uniqueRatsAtSuspiciousAge);

%% Session-level Summary for reversal_trial (for reference)
reversalTrials = suspiciousData.reversal_trial;
fprintf('\nSession-level Reversal Trial Summary:\n');
fprintf('  Min: %.2f\n', min(reversalTrials));
fprintf('  Max: %.2f\n', max(reversalTrials));
fprintf('  Mean: %.2f\n', mean(reversalTrials, 'omitnan'));
fprintf('  Std: %.2f\n', std(reversalTrials, 'omitnan'));

%% Within-Rat Averaging: Compute rat-level average reversal trial number
[G_rats, theseRats] = findgroups(suspiciousData.rat);
ratMeanReversal = splitapply(@(x) mean(x, 'omitnan'), suspiciousData.reversal_trial, G_rats);

overallRatMean = mean(ratMeanReversal, 'omitnan');
fprintf('\nOverall rat-level mean reversal trial (within-rat averaging): %.2f\n', overallRatMean);

disp('Rat-level reversal trial numbers:');
T_rats = table(theseRats, ratMeanReversal, 'VariableNames', {'Rat','MeanReversalTrial'});
disp(T_rats);

%% Visualize the distribution of rat-level reversal trial numbers
figure('Name','Investigation at Age 120: Rat-level Averages');
subplot(1,2,1);
boxplot(ratMeanReversal);
xlabel('All Rats');
ylabel('Mean Reversal Trial Number');
title('Boxplot of Rat-level Mean Reversal Trial Number');
grid on;

subplot(1,2,2);
histogram(ratMeanReversal, 'BinWidth', 10);
xlabel('Mean Reversal Trial Number');
ylabel('Count');
title('Histogram of Rat-level Mean Reversal Trial Number');
grid on;