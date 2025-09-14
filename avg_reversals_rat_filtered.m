% Input/output files
inputFile = '/Users/olly/Desktop/ASI/ASI/MAT Files/reversals_per_session.csv';
outputFile = '/Users/olly/Desktop/ASI/ASI/MAT Files/avg_reversals_per_rat.csv';
coreDataFile = '/Users/olly/Desktop/ASI/ASI/MAT Files/core_ASI_data.csv';

% Load session-level data and core metadata
data = readtable(inputFile);
core = readtable(coreDataFile);

% Get unique rats and ages
unique_rats = unique(data.Rat);
unique_ages = sort(unique(data.Age));

% Create column names
age_colnames = "AvgReversals_P" + string(unique_ages);
var_names = ["Rat"; "Group"; age_colnames(:)];

% Create result table
nRats = numel(unique_rats);
nAges = numel(unique_ages);
result = table('Size', [nRats, 2 + nAges], ...
               'VariableTypes', ["string", "string", repmat("double", 1, nAges)], ...
               'VariableNames', cellstr(var_names));

% Loop over rats and calculate average reversals per age (after filtering)
for r = 1:nRats
    rat_id = unique_rats{r};
    result.Rat(r) = rat_id;

    % Lookup group info from core metadata
    match_idx = strcmp(core.Rat, rat_id);
    if any(match_idx)
        result.Group(r) = core.Group(find(match_idx, 1));  % pick first match
    else
        result.Group(r) = "UNKNOWN";
        warning("Rat %s not found in core_ASI_data.csv", rat_id);
    end

    for a = 1:nAges
        age_val = unique_ages(a);

        % Get all session reversal counts for this rat at this age
        mask = strcmp(data.Rat, rat_id) & data.Age == age_val;
        revs = data.Reversals(mask);

        if isempty(revs)
            continue;
        end

        % IQR filtering
        Q1 = quantile(revs, 0.25);
        Q3 = quantile(revs, 0.75);
        IQR_val = Q3 - Q1;
        lower = Q1 - 1.5 * IQR_val;
        upper = Q3 + 1.5 * IQR_val;
        revs_filtered = revs(revs >= lower & revs <= upper);

        if ~isempty(revs_filtered)
            colname = "AvgReversals_P" + string(age_val);
            result{r, colname} = mean(revs_filtered);
        end
    end
end

% Save final result table
writetable(result, outputFile);
fprintf('✅ Saved per-rat reversal averages + group to:\n%s\n', outputFile);