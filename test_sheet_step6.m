%% export_RPS_metrics_singleRowPerRat.m
% This script reads core_ASI_data.csv and computes the average reversals per
% session for each rat (averaging across sessions). It uses the columns:
% "Reversals_per_Session_Age_<age>" for each age defined in the 'ages' vector.
% The output is a spreadsheet with one row per rat and one column per age.
% The columns are named RPS_Age_<age> (e.g., RPS_Age_30).
%
% Adjustable parameters are defined at the top.

%% === Adjustable Parameters ===
csvFile = '/Users/olly/Desktop/ASI/ASI/MAT Files/core_ASI_data.csv';

% Age bins and corresponding column names in the data
ages = [30, 50, 70, 90, 120, 190];
revPerSessCols = strcat("Reversals_per_Session_Age_", string(ages));

% Group filtering: specify which group(s) to include.
groupList = {'SOC'};  % set to [] to include all groups

%% === Load Data ===
T = readtable(csvFile);

% Filter by group if specified
if ~isempty(groupList)
    T = T(ismember(T.Group, groupList), :);
end

%% === Get Unique Rats ===
% Assumes the rat identifier is stored in the column 'Rat'
uniqueRats = unique(T.Rat, 'stable');
nRats = numel(uniqueRats);
nAges = numel(ages);

%% === Compute Average Reversals per Session for Each Rat ===
% Preallocate matrix for averaged values
RPS = NaN(nRats, nAges);

for i = 1:nRats
    % Get all sessions for the current rat
    ratData = T(strcmp(T.Rat, uniqueRats{i}), :);
    for j = 1:nAges
        colName = revPerSessCols{j};  % e.g., "Reversals_per_Session_Age_30"
        if ismember(colName, T.Properties.VariableNames)
            values = ratData.(colName);
            % Remove NaNs and compute the mean
            values = values(~isnan(values));
            if ~isempty(values)
                RPS(i, j) = mean(values);
            end
        else
            warning("Missing column: %s", colName);
        end
    end
end

%% === Build the Output Table ===
% Create column names for the output: e.g., RPS_Age_30, RPS_Age_50, ...
outColNames = arrayfun(@(a) sprintf('RPS_P%d', a), ages, 'UniformOutput', false);

% Create a table with one row per rat
T_out = table(uniqueRats, 'VariableNames', {'rat'});
T_out = [T_out, array2table(RPS, 'VariableNames', outColNames)];

%% === Write the Output Table to CSV ===
outputFile = '/Users/olly/Desktop/ASI/ASI/MAT Files/RPS_test.csv';
writetable(T_out, outputFile);
fprintf('Exported %d rows (one per rat) to %s\n', height(T_out), outputFile);