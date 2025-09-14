
processDir = '/Users/oliver/Desktop/ASI/ASI/MAT Files/raw_cpdata_190'
age = 190

process1_cp_data(processDir, age)


function process1_cp_data(directory, age)
    % Processes CP data by adding delta x, delta y, and slope columns.
    % Searches for all CP CSV files in the given age directory.

    % Use the given directory directly
    cpDataDir = directory;
    if ~isfolder(cpDataDir)
        disp(['Directory not found: ', cpDataDir]);
        return;
    end

    % Define the file pattern for CP data files
    filePattern = fullfile(cpDataDir, 'CP_*.csv');
    files = dir(filePattern);

    % Check if any files were found
    if isempty(files)
        disp('No CP CSV files found.');
        return;
    end

    % Loop through each CP CSV file
    for n = 1:length(files)
        csvFile = fullfile(files(n).folder, files(n).name);
        disp(['Processing file: ', files(n).name]);

        % Read the CP data
        cpData = readtable(csvFile);

        % Compute delta x (change in trial number), first row vs. 0
        deltaX = [cpData.Trial(1); diff(cpData.Trial)];

        % Compute delta y (change in cumulative value), first row vs. 0
        deltaY = [cpData.Cumulative_Value(1); diff(cpData.Cumulative_Value)];

        % Compute slope (deltaY / deltaX) safely (avoid division by zero)
        slope = deltaY ./ deltaX;
        slope(deltaX == 0) = 0; % Avoid NaN due to division by zero

        % Append new columns to the table
        cpData.Delta_X = deltaX;
        cpData.Delta_Y = deltaY;
        cpData.Slope = slope;

        % Save the updated table back to the CSV file
        writetable(cpData, csvFile);
        disp(['Updated file saved: ', csvFile]);
    end

    disp('Finished processing all CP files.');
end
