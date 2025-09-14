% File paths
exploreExploitFile = '/Users/oliver/Desktop/ASI/ASI/MAT Files/exploit_explore_data.csv';

% Read data
data = readtable(exploreExploitFile);

% Compute "CP1 After - Reversal Trial" (First CP latency)
data.CP1_Latency = data.CP1_Trial_After - data.Reversal_Trial;

% Identify cases where no CP was detected after reversal (NaNs)
data.No_CP_After = isnan(data.CP1_Latency);

% Define groups and colors
groups = unique(data.Group);
colors = lines(length(groups)); % Generate distinct colors for groups

% Get unique ages
ages = unique(data.Age);

% Compute proportion of NaNs (cases where no CP was detected)
nanProportions = arrayfun(@(age) mean(data.No_CP_After(data.Age == age)), ages);

% Initialize figure
figure;
yyaxis left;
hold on;

% Plot mean latency per age for each group
for g = 1:length(groups)
    groupIdx = strcmp(data.Group, groups{g});
    meanLatencies = arrayfun(@(age) mean(data.CP1_Latency(groupIdx & data.Age == age), 'omitnan'), ages);
    plot(ages, meanLatencies, '-o', 'Color', colors(g, :), 'LineWidth', 2, 'MarkerFaceColor', colors(g, :), 'DisplayName', groups{g});
end

ylabel('Mean Trials to First CP After Reversal');
ylim([0, max(data.CP1_Latency, [], 'omitnan') + 5]); % Adjust Y-axis range

% Overlay proportion of NaNs as a bar plot on the right axis
yyaxis right;
bar(ages, nanProportions, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
ylabel('Proportion of Reversals with No CP After');

% Labels and title
xlabel('Age');
title('Mean Latency to First Change Point After Reversal & Missing CPs');
legend('Location', 'best');

grid on;
hold off;

disp('Plot generated: Mean Latency to First CP After Reversal with Proportion of NaNs');