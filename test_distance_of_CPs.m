% Parameters
critVal = 3.8; % Set your crit value
threshold = 2; % Changeable threshold for minimum trial distance
dataDir = '/Users/oliver/Desktop/ASI/ASI/MAT Files/';
dataFile = sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal);
filePath = fullfile(dataDir, dataFile);

% Load data
data = readtable(filePath);
data.Properties.VariableNames = lower(data.Properties.VariableNames);

% Initialize counter
totalCloseCPs = 0;
closeCPsByAge = zeros(length(unique(data.age)), 1);
ages = unique(data.age);

% Loop through each session
for i = 1:height(data)
    % Collect CP trials before and after
    cpTrials = [data.cp1_trial_before(i), data.cp2_trial_before(i), data.cp3_trial_before(i), data.cp4_trial_before(i), ...
                data.cp1_trial_after(i), data.cp2_trial_after(i), data.cp3_trial_after(i), data.cp4_trial_after(i)];
    cpTrials = cpTrials(~isnan(cpTrials)); % Remove NaNs
    
    % Calculate differences
    trialDiffs = diff(cpTrials);
    
    % Count close CPs
    closeCount = sum(trialDiffs < threshold);
    totalCloseCPs = totalCloseCPs + closeCount;
    
    % Count per age
    ageIdx = find(ages == data.age(i));
    closeCPsByAge(ageIdx) = closeCPsByAge(ageIdx) + closeCount;
end

% Display result
disp(['Total number of CPs fewer than ', num2str(threshold), ' trials apart: ', num2str(totalCloseCPs)]);

% Display by age
for i = 1:length(ages)
    disp(['Age ', num2str(ages(i)), ': ', num2str(closeCPsByAge(i)), ' close CPs']);
end