condensed_segment_cp_extract(3, true)


function condensed_segment_cp_extract(critVal, showPlots)

% ----------------- PARAMETERS --------------------
% Base directories
sharedSegmentDir = '/Users/olly/Desktop/ASI/ASI/MAT Files/segment4.15';  % All segments and REV markers
baseOutputDir = '/Users/olly/Desktop/ASI/ASI/MAT Files/new_raw_CP_4.15';  % <- NEW SAVE FOLDER

% Ages to process
ages = [190];

% Create output folder if it doesn't exist
if ~exist(baseOutputDir, 'dir')
    mkdir(baseOutputDir);
end

% Output folder name based on Crit value
cpDataDir = fullfile(baseOutputDir, sprintf('raw_cpdata_crit%.1f', critVal));
if ~exist(cpDataDir, 'dir')
    mkdir(cpDataDir);
end
disp(['CP Data will be saved to: ', cpDataDir]);
% --------------------------------------------------

for a = 1:length(ages)
    age = ages(a);
    disp(['Processing age: ', num2str(age)]);
    
    % Use same folder for both segments and reversal markers
    segmentDir = sharedSegmentDir;
    
    % Locate REV_MARKERS file
    markerFile = fullfile(sharedSegmentDir, sprintf('REV_MARKERS_P%d.txt', age));
    if ~isfile(markerFile)
        disp(['Reversal marker file not found: ', markerFile]);
        continue;
    end
    markerData = readtable(markerFile, 'Delimiter', 'tab');
    
    % Find all .txt segment files for this age
    filePattern = fullfile(sharedSegmentDir, sprintf('SI_DEV*_P%d_PROB_*_REV*_REV*.txt', age));
    files = dir(filePattern);
    
    if isempty(files)
        disp(['No segmented .TXT files found for age ', num2str(age)]);
        continue;
    end
    
    % Loop through each file
    for n = 1:length(files)
        txtFile = fullfile(files(n).folder, files(n).name);
        disp(['Processing file: ' files(n).name]);

        % Load binary data
        binaryData = readmatrix(txtFile);
        binaryData = double(binaryData);
        if isempty(binaryData)
            disp(['Skipping empty file: ' files(n).name]);
            continue;
        end
        
        % Run cp_wrapper
        cp_wrapper(binaryData, 1, 4, critVal);  % Test 4 = Chi-square
        
        % Extract CP data from figure
        figs = findall(0, 'Type', 'Figure'); 
        latestFig = max([figs.Number]);
        figure(latestFig);

        ax3 = subplot(3,1,3); 
        textObjs = findall(ax3, 'Type', 'Text'); 
        textStr = textObjs.String; 

        lines = split(textStr, '\n'); 
        cpData = [];
        for i = 2:length(lines) 
            tokens = regexp(lines{i}, 'CP (\d+): \(Trial (\d+), Cumulative ([\d.]+)\)', 'tokens');
            if ~isempty(tokens)
                cpRow = str2double(tokens{1});
                cpData = [cpData; cpRow]; 
            end
        end

        % Save CP data
        [~, baseName, ~] = fileparts(files(n).name);
        outputFile = fullfile(cpDataDir, sprintf('CP_%s.csv', baseName));
        cp_table = array2table(cpData, 'VariableNames', {'CP_Number', 'Trial', 'Cumulative_Value'});
        writetable(cp_table, outputFile);
        disp(['Saved CP data to ', outputFile]);
        
        % Add Reversal Marker
        ax1 = subplot(3,1,1); 
        hold on;
        rowIdx = find(strcmp(markerData.Segment_File, files(n).name));
        if ~isempty(rowIdx)
            relativeReversal = markerData.Relative_Middle_Reversal(rowIdx);
            xline(relativeReversal, 'r--', 'LineWidth', 1.5);
            text(relativeReversal, max(ylim)-0.1*(max(ylim)-min(ylim)), 'Reversal', ...
                'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
            disp(['Added reversal marker at trial ', num2str(relativeReversal)]);
        end
        hold off;
        
        % Optionally close plot
        if showPlots == false
            close(latestFig);
        end
    end
end

disp('✅ Finished processing all files.');

end