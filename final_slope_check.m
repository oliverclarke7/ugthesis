% File path
exploreExploitFile = '/Users/oliver/Desktop/ASI/ASI/MAT Files/exploit_explore_data.csv';

% Load data
data = readtable(exploreExploitFile);

% Initialize counter
countCloseCPs = 0;

% Loop through each row
for i = 1:height(data)
    % Extract CP trial numbers (ignore NaNs)
    CPs_before = [data.CP1_Trial_Before(i), data.CP2_Trial_Before(i), data.CP3_Trial_Before(i)];
    CPs_before = CPs_before(~isnan(CPs_before));
    
    CPs_after = [data.CP1_Trial_After(i), data.CP2_Trial_After(i), data.CP3_Trial_After(i)];
    CPs_after = CPs_after(~isnan(CPs_after));
    
    % Check before reversal CPs
    if length(CPs_before) >= 2
        diffs = diff(CPs_before);
        if any(diffs < 2)
            countCloseCPs = countCloseCPs + 1;
            fprintf('Row %d: CPs BEFORE reversal less than 3 trials apart\n', i);
        end
    end
    
    % Check after reversal CPs
    if length(CPs_after) >= 2
        diffs = diff(CPs_after);
        if any(diffs < 2)
            countCloseCPs = countCloseCPs + 1;
            fprintf('Row %d: CPs AFTER reversal less than 3 trials apart\n', i);
        end
    end
end

% Final count
fprintf('\nTotal cases with CPs less than 3 trials apart: %d\n', countCloseCPs);