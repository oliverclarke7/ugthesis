% === File Path ===
csvFile = '/Users/olly/Desktop/ASI/ASI/MAT Files/core_ASI_data.csv';

% === Age bins and corresponding reversal/session columns ===
ages = [30, 50, 70, 90, 120, 190];
revPerSessCols = strcat("Reversals_per_Session_Age_", string(ages));
groupList = {'SOC'};

% === Load Data ===
T = readtable(csvFile);

% === Initialize matrices for mean and SEM ===
meanRevPerSess = zeros(length(groupList), length(ages));
semRevPerSess  = zeros(length(groupList), length(ages));

% === Compute mean and SEM for each group/age ===
for g = 1:length(groupList)
    groupName = groupList{g};
    groupData = T(strcmp(T.Group, groupName), :);
    
    for a = 1:length(ages)
        colName = revPerSessCols{a};
        if ismember(colName, T.Properties.VariableNames)
            values = groupData.(colName);
            values = values(~isnan(values));  % Remove NaNs
            n = numel(values);
            meanRevPerSess(g, a) = mean(values);
            semRevPerSess(g, a)  = std(values) / sqrt(n);
        else
            warning("Missing column: %s", colName);
        end
    end
end

% === Plotting ===
figure; hold on;
colors = lines(length(groupList)); % Distinct colors for each group

for g = 1:length(groupList)
    % Plot line with error bars
    errorbar(ages, meanRevPerSess(g,:), semRevPerSess(g,:), ...
        '-o', 'LineWidth', 2, 'DisplayName', groupList{g}, 'Color', colors(g,:));
end

xlabel('Age (Days)');
ylabel('Average Reversals per Session');
title('Average Reversals per Session by Age and Group (with SEM)');
legend('Location', 'northwest');
grid on;