% Set path to PROB directory
probFolder = '/Users/olly/Desktop/ASI/ASI/MAT Files/PROB';
saveFolder = '/Users/olly/Desktop/ASI/ASI/MAT Files'
outputFile = fullfile(saveFolder, 'reversals_per_session.csv');

files = dir(fullfile(probFolder, '*.MAT'));

% Initialize output storage
ratNames = {};
ages = [];
dates = {};
reversalCounts = [];

reversal_values = [14, 6, 2];  % valid values only

for f = 1:length(files)
    filename = files(f).name;
    fullpath = fullfile(probFolder, filename);

    % Extract metadata from filename
    tokens = regexp(filename, 'SI_(DEV\d+)_P(\d+)_PROB_(\d+)\.MAT', 'tokens', 'ignorecase');
    if isempty(tokens)
        warning("Skipped file (bad name): %s", filename);
        continue;
    end
    tokens = tokens{1};
    rat = tokens{1};
    age = str2double(tokens{2});
    date = tokens{3};

    % Load .MAT session file
    raw = load(fullpath);
    fn = fieldnames(raw);
    data = raw.(fn{1});

    if size(data,1) < 5
        warning("File has <5 rows, skipping: %s", filename);
        continue;
    end

    % Extract rows 3–5 (port probabilities)
    probs = data(3:5, 4:end);  % start from column 4

    % Count valid state transitions
    rev_count = 0;
    for t = 2:size(probs, 2)
        prev = probs(:, t-1);
        curr = probs(:, t);

        % Only count if all values are in [14, 6, 2] (i.e., valid state)
        if all(ismember(prev, reversal_values)) && all(ismember(curr, reversal_values))
            if ~isequal(prev, curr)
                rev_count = rev_count + 1;
            end
        end
    end

    % Store results
    ratNames{end+1,1} = rat;
    ages(end+1,1) = age;
    dates{end+1,1} = date;
    reversalCounts(end+1,1) = rev_count;
end

%% --- Outlier Detection using IQR ---

Q1 = quantile(reversalCounts, 0.25);
Q3 = quantile(reversalCounts, 0.75);
IQR_val = Q3 - Q1;
lower_bound = Q1 - 1.5 * IQR_val;
upper_bound = Q3 + 1.5 * IQR_val;
isOutlier = (reversalCounts < lower_bound) | (reversalCounts > upper_bound);

%% --- Save to CSV ---

T = table(ratNames, ages, dates, reversalCounts, isOutlier, ...
    'VariableNames', {'Rat', 'Age', 'Date', 'Reversals', 'Outlier'});

writetable(T, outputFile);
fprintf('✅ Cleaned reversal data (corrected logic) saved to:\n%s\n', outputFile);