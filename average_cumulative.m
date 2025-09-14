plot_fixed_window_cumulative_curve(190)

function plot_fixed_window_cumulative_curve(age)
    % 0) USER PARAMETERS
    %    Adjust directory paths, file search pattern, etc. as needed
    baseSegmentDir = '/Users/olly/Desktop/ASI/ASI/MAT Files/segmented data';
    baseMarkerDir  = '/Users/olly/Desktop/ASI/ASI/MAT Files/segmented data';
    
    % Window size around reversal
    preTrials  = 50;  % how many trials before the reversal to plot
    postTrials = 50;  % how many trials after the reversal to plot
    
    % The x-axis from -100 up to +100 in steps of 1
    xWindow = -preTrials:postTrials;
    xLen    = length(xWindow);  % should be 201 if pre=100 and post=100
    
    % 1) Locate the folder for the given age
    segmentDir = fullfile(baseSegmentDir, sprintf('segment_%d', age));
    if ~exist(segmentDir, 'dir')
        error('Segment folder not found for age %d', age);
    end
    
    % 2) Load the reversal markers
    markerFile = fullfile(baseMarkerDir, sprintf('REV_MARKERS_P%d.txt', age));
    if ~isfile(markerFile)
        error('Reversal marker file not found: %s', markerFile);
    end
    markerData = readtable(markerFile, 'Delimiter', 'tab');  %#ok<*DTXTRD>
    
    % 3) Find all relevant .txt files
    filePattern = fullfile(segmentDir, sprintf('SI_DEV*_P%d_PROB_*_REV*_REV*.txt', age));
    files = dir(filePattern);
    if isempty(files)
        error('No segmented .TXT files found for age %d', age);
    end
    
    % Initialize storage for all rats' curves
    allCurves = [];  % each row => one file, each column => cumsum at window index
    
    % 4) Process each file
    for f = 1:length(files)
        fname = files(f).name;
        fpath = fullfile(files(f).folder, fname);
        disp(['Processing file: ' fname]);
        
        % 4a) Load binary data
        binaryData = readmatrix(fpath);
        if isempty(binaryData)
            warning('Empty file: %s, skipping...', fname);
            continue;
        end
        binaryData = double(binaryData);
        T = length(binaryData);
        
        % 4b) Find the reversal trial R for this file
        %     We assume 'Segment_File' in markerData matches the filename.
        rowIdx = find(strcmp(markerData.Segment_File, fname));
        if isempty(rowIdx)
            warning('No reversal marker found for file: %s, skipping...', fname);
            continue;
        end
        reversalTrial = markerData.Relative_Middle_Reversal(rowIdx);
        
        % 4c) Compute the cumsum
        csumVals = cumsum(binaryData);
        
        % 4d) Create a per-file curve aligned around R
        %     xWindow(i) => relative offset. We want csumVals(R + xWindow(i)).
        fileCurve = nan(1, xLen);
        for i = 1:xLen
            trialIndex = reversalTrial + xWindow(i); 
            if trialIndex >= 1 && trialIndex <= T
                fileCurve(i) = csumVals(trialIndex);
            end
        end
        
        % Store in allCurves
        allCurves = [allCurves; fileCurve];
    end
    
    % 5) Compute average curve (ignore NaNs)
    avgCurve = nanmean(allCurves, 1);
    
    % 6) Plot
    figure('Color','w','Name',sprintf('Age %d', age));
    hold on; box on;
    
    % Light gray lines for individual curves
    for r = 1:size(allCurves,1)
        plot(xWindow, allCurves(r,:), 'Color', [0.8 0.8 0.8], 'LineStyle', '-');
    end
    
    % Bold black line for the average
    plot(xWindow, avgCurve, 'k', 'LineWidth', 2.5);
    
    % Vertical red dashed line at x=0 => Reversal
    xline(0, 'r--', 'LineWidth', 1.5);
    
    title(sprintf('Averaged Cumulative Curve (Age = %d)', age), 'FontSize', 12);
    xlabel(sprintf('Trials Relative to Reversal (0 = Reversal)'));
    ylabel('Average Cumulative Value');
    axis tight;
    grid on;
    
    disp('✅ Finished plotting cumulative curves');
end