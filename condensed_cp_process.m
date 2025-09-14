

condensed_cp_process_allcrit('/Users/olly/Desktop/ASI/ASI/MAT Files/new_raw_CP_4.15', 1.5)

function condensed_cp_process_allcrit(baseCritDir, critVal)
% Processes all CP CSV files across all ages for a specific Crit value.
% Adds Delta_X, Delta_Y, and Slope columns.

% Define ages to process
ages = [30, 50, 70, 90, 120, 190];

for a = 1:length(ages)
    age = ages(a);
    disp(['Processing CP files for age: ', num2str(age)]);
    
    % CP data directory based on critVal
    cpDataDir = fullfile(baseCritDir, sprintf('raw_cpdata_crit%.1f', critVal));
    
    if ~isfolder(cpDataDir)
        disp(['Directory not found: ', cpDataDir]);
        continue;
    end

    % Locate CP CSVs
    filePattern = fullfile(cpDataDir, sprintf('CP_SI_DEV*_P%d_PROB_*.csv', age));
    files = dir(filePattern);
    
    if isempty(files)
        disp(['No CP CSV files found for age ', num2str(age)]);
        continue;
    end

    % Loop through each CP CSV
    for n = 1:length(files)
        csvFile = fullfile(files(n).folder, files(n).name);
        disp(['Processing file: ' files(n).name]);

        cpData = readtable(csvFile);

        % Calculate delta values
        deltaX = [cpData.Trial(1); diff(cpData.Trial)];
        deltaY = [cpData.Cumulative_Value(1); diff(cpData.Cumulative_Value)];
        slope = deltaY ./ deltaX;
        slope(deltaX == 0) = 0;  % Avoid division by zero

        % Append to table
        cpData.Delta_X = deltaX;
        cpData.Delta_Y = deltaY;
        cpData.Slope = slope;

        % Save updated file
        writetable(cpData, csvFile);
        disp(['Updated file saved: ', csvFile]);
    end
end

disp('✅ Finished processing all CP files across all ages.');

end