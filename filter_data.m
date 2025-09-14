baseCritDir = '/Users/olly/Desktop/ASI/ASI/MAT Files/';
critVal = 3.2;

filter_cp_remove_dummy(baseCritDir, critVal);

function filter_cp_remove_dummy(baseCritDir, critVal)
% Remove dummy CPs (last CP with NaN slope) in exploit_explore_data_critX.X.csv
% Does NOT shift columns, keeps structure intact
% Outputs _filtered CSV

inputFile = fullfile(baseCritDir, sprintf('exploit_explore_data_crit%.1f.csv', critVal));
outputFile = fullfile(baseCritDir, sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal));

% Load
T = readtable(inputFile);
T.Properties.VariableNames = lower(T.Properties.VariableNames);

for i = 1:height(T)
    % ----- Before -----
    for p = 1:4
        if isnan(T.(sprintf('cp%d_slope_before', p))(i))
            T.(sprintf('cp%d_trial_before', p))(i) = NaN;
            T.(sprintf('cp%d_cumulative_before', p))(i) = NaN;
        end
    end
    
    % ----- After -----
    for p = 1:4
        if isnan(T.(sprintf('cp%d_slope_after', p))(i))
            T.(sprintf('cp%d_trial_after', p))(i) = NaN;
            T.(sprintf('cp%d_cumulative_after', p))(i) = NaN;
        end
    end
end

% Save
writetable(T, outputFile);
disp(['✅ Dummy CP rows removed and saved to: ', outputFile]);

end