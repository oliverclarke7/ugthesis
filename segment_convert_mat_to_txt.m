function segment_convert_mat_to_txt(directory, savedir, age)
    % Converts .MAT files to .TXT format for each valid reversal segment.
    % Stores all reversal markers in a single file: REV_MARKERS_P[age].txt.
    
    % Define the file pattern
    filePattern = ['SI_DEV*_P' num2str(age) '_PROB_*.mat'];
    d = dir(fullfile(directory, filePattern));
    
    % Ensure the output directory exists
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end

    % Initialize storage for all reversal markers
    reversal_metadata = {}; 

    for n = 1:length(d)
        sname = d(n).name;
        filePath = fullfile(directory, sname);
        
        % Load the .MAT file
        load(filePath, 'data');  

        % Extract choices and reward probabilities
        choices = data(1, 4:end);  % Row 1 (Trial choices)
        probabilities = data(3:5, 4:end);  % Rows 3-5 (Reward probabilities)

        % Identify reversal points (columns where probabilities change)
        reversals = find(any(diff(probabilities, 1, 2) ~= 0, 1)) + 4; % Adjust for MATLAB indexing

        % Ensure there are enough reversals for segmentation
        if length(reversals) < 2
            fprintf('Skipping %s: Not enough reversals.\n', sname);
            continue;
        end

        % First file (special case: trial 1 → Reversal 2)
        startIdx = 1;
        endIdx = reversals(2) - 1;
        middleReversal = reversals(1);
        relativeReversal = middleReversal - startIdx + 1;

        [~, maxIdx] = max(probabilities(:, middleReversal - 4));
        bestNP = maxIdx + 1; % Convert row index to NP (since rows map to NP2, NP3, NP4)
        segmentChoices = choices(startIdx:endIdx);
        correct_choices = (segmentChoices == bestNP); % 1 if best NP chosen, 0 otherwise

        [~, baseName, ~] = fileparts(sname);
        txtFileName = sprintf('%s_REV1_REV2.txt', baseName);
        writematrix(correct_choices', fullfile(savedir, txtFileName), 'Delimiter', 'tab');
        fprintf('Converted segment 1-2 of %s to %s\n', sname, txtFileName);

        % Store metadata for first segment
        reversal_metadata{end+1, 1} = txtFileName;
        reversal_metadata{end, 2} = relativeReversal;

        % Iterate through remaining valid reversal segments
        for i = 2:length(reversals) - 1
            startIdx = reversals(i - 1);
            endIdx = reversals(i + 1) - 1;
            middleReversal = reversals(i);
            relativeReversal = middleReversal - startIdx + 1;

            % Identify NP with the highest reward probability BEFORE the middle reversal
            [~, maxIdx] = max(probabilities(:, middleReversal - 4));
            bestNP = maxIdx + 1; % Convert row index to NP (since rows map to NP2, NP3, NP4)

            % Extract the choice data for this segment
            segmentChoices = choices(startIdx:endIdx);
            correct_choices = (segmentChoices == bestNP); % 1 if best NP chosen, 0 otherwise

            % Create the output filename
            txtFileName = sprintf('%s_REV%d_REV%d.txt', baseName, i - 1, i + 1);
            writematrix(correct_choices', fullfile(savedir, txtFileName), 'Delimiter', 'tab');

            fprintf('Converted segment %d-%d of %s to %s\n', i - 1, i + 1, sname, txtFileName);

            % Store metadata
            reversal_metadata{end+1, 1} = txtFileName;
            reversal_metadata{end, 2} = relativeReversal;
        end
    end

    % Save reversal metadata to a single file for all segments of this age
    metaFileName = fullfile(savedir, sprintf('REV_MARKERS_P%d.txt', age));
    writetable(cell2table(reversal_metadata, 'VariableNames', {'Segment_File', 'Relative_Middle_Reversal'}), ...
        metaFileName, 'Delimiter', 'tab');

    fprintf('Saved all reversal markers for age %d to %s\n', age, metaFileName);
end