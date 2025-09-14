coreDataPath = '/Users/oliver/Desktop/ASI/ASI/MAT Files/core_ASI_data.csv';
dataDirectory = '/Users/oliver/Desktop/ASI/ASI/MAT Files';
probDirectory = '/Users/oliver/Desktop/ASI/ASI/MAT Files/PROB';
ages = [30, 50, 70, 90, 120, 190]; % Add all relevant ages

append_reversal_counts(coreDataPath, dataDirectory, probDirectory, ages);

function append_reversal_counts(coreDataPath, dataDirectory, probDirectory, ages)
    % Appends reversal counts and reversals per session for each age to the core ASI data.
    % Counts how many segment files exist per rat in each raw_cpdata_age folder.
    % Also counts the number of sessions per rat by checking the number of PROB files.

    % Read the core ASI data
    coreData = readtable(coreDataPath);
    
    % Iterate over each specified age
    for i = 1:length(ages)
        age = ages(i);
        
        % Define the directory for this age's CP data
        cpDataDir = fullfile(dataDirectory, sprintf('raw_cpdata_%d', age));
        if ~isfolder(cpDataDir)
            disp(['Directory not found: ', cpDataDir]);
            continue;
        end
        
        % Define the directory for PROB files
        probDataDir = fullfile(probDirectory);
        if ~isfolder(probDataDir)
            disp(['Directory not found: ', probDataDir]);
            continue;
        end
        
        % Get list of all CP files
        filePattern = fullfile(cpDataDir, 'CP_*.csv');
        files = dir(filePattern);
        
        % Extract rat names from filenames, removing "SI_"
        ratIDs = regexp({files.name}, 'CP_SI_(DEV\d+)_P', 'tokens', 'once');
        ratIDs = [ratIDs{:}]; % Convert from cell array of cells to flat cell array
        
        % Count occurrences of each rat ID
        uniqueRats = unique(ratIDs);
        reversalCounts = zeros(height(coreData), 1);
        sessionCounts = zeros(height(coreData), 1);
        
        for j = 1:length(uniqueRats)
            ratName = uniqueRats{j};
            count = sum(strcmp(ratIDs, ratName));
            
            % Find matching row in core data
            rowIdx = find(strcmp(coreData.Rat, ratName));
            if ~isempty(rowIdx)
                reversalCounts(rowIdx) = count;
            end
            
            % Count the number of session files for this rat
            sessionPattern = fullfile(probDataDir, sprintf('SI_%s_P%d_PROB_*.mat', ratName, age));
            sessionFiles = dir(sessionPattern);
            sessionCounts(rowIdx) = length(sessionFiles);
        end
        
        % Compute reversals per session
        reversalsPerSession = reversalCounts ./ sessionCounts;
        reversalsPerSession(isnan(reversalsPerSession)) = 0; % Handle division by zero
        
        % Add new columns to core data
        columnReversal = sprintf('Reversals_Age_%d', age);
        columnReversalSession = sprintf('Reversals_per_Session_Age_%d', age);
        coreData.(columnReversal) = reversalCounts;
        coreData.(columnReversalSession) = reversalsPerSession;
    end
    
    % Save updated core data
    writetable(coreData, coreDataPath);
    disp('Updated core ASI data saved.');
end