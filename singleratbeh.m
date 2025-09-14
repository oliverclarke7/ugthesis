condensed_segment_cp_extract_rat(3.2, true, 'DEV021016007');

function condensed_segment_cp_extract_rat(critVal, showPlots, ratID)
    % ----------------- PARAMETERS --------------------
    % Base directories
    baseSegmentDir = '/Users/olly/Desktop/ASI/ASI/MAT Files/segmented data';
    baseOutputDir = '/Users/olly/Desktop/ASI/ASI/MAT Files/';

    % Ages to process (can be modified as needed)
        ages = [190];

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

        segmentDir = fullfile(baseSegmentDir, sprintf('segment_%d', age));

        % Check if folder exists
        if ~exist(segmentDir, 'dir')
            disp(['Segment folder not found for age ', num2str(age)]);
            continue;
        end

        % Locate REV_MARKERS file
        markerFile = fullfile(baseSegmentDir, sprintf('REV_MARKERS_P%d.txt', age));
        if ~isfile(markerFile)
            disp(['Reversal marker file not found: ', markerFile]);
            continue;
        end
        markerData = readtable(markerFile, 'Delimiter', 'tab');

        % Find segmented .txt files with given age
        filePattern = fullfile(segmentDir, sprintf('SI_DEV*_P%d_PROB_*_REV*_REV*.txt', age));
        files = dir(filePattern);

        % Filter files for the specified rat if ratID is provided
        if nargin >= 3 && ~isempty(ratID)
            files = files(contains({files.name}, ratID));
        end

        if isempty(files)
            disp(['No segmented .TXT files found for age ', num2str(age), ' and rat ', ratID]);
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

            % Run cp_wrapper with new Crit value
            cp_wrapper(binaryData, 1, 4, critVal);  % Test 4 = Chi-square, binary data

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

            % Save CP Data
            [~, baseName, ~] = fileparts(files(n).name);
            outputFile = fullfile(cpDataDir, sprintf('CP_%s.csv', baseName));
            cp_table = array2table(cpData, 'VariableNames', {'CP_Number', 'Trial', 'Cumulative_Value'});
            writetable(cp_table, outputFile);
            disp(['Saved CP data to ', outputFile]);

            % Add Reversal Marker to the first subplot
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