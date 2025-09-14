coreDataPath = '/Users/oliver/Desktop/ASI/ASI/MAT Files/core_ASI_data.csv';
ages = [30, 50, 70, 90, 120, 190]; % Modify as needed

plot_avg_reversals(coreDataPath, ages);

function plot_avg_reversals(coreDataPath, ages)
    % Plots average reversals and reversals per session by age.

    % Read the core ASI data
    coreData = readtable(coreDataPath);
    
    % Initialize storage for averages
    avgReversals = zeros(length(ages), 1);
    avgReversalsPerSession = zeros(length(ages), 1);
    
    % Loop through each age to compute means
    for i = 1:length(ages)
        age = ages(i);
        columnReversal = sprintf('Reversals_Age_%d', age);
        columnReversalSession = sprintf('Reversals_per_Session_Age_%d', age);
        
        if ismember(columnReversal, coreData.Properties.VariableNames)
            avgReversals(i) = mean(coreData.(columnReversal), 'omitnan');
        else
            avgReversals(i) = NaN;
        end
        
        if ismember(columnReversalSession, coreData.Properties.VariableNames)
            avgReversalsPerSession(i) = mean(coreData.(columnReversalSession), 'omitnan');
        else
            avgReversalsPerSession(i) = NaN;
        end
    end
    
    % Plot average reversals by age
    figure;
    subplot(1,2,1);
    bar(ages, avgReversals, 'FaceColor', [0.2, 0.6, 0.8]);
    xlabel('Age');
    ylabel('Average Reversals');
    title('Average Reversals by Age');
    grid on;
    
    % Plot average reversals per session by age
    subplot(1,2,2);
    bar(ages, avgReversalsPerSession, 'FaceColor', [0.8, 0.4, 0.2]);
    xlabel('Age');
    ylabel('Average Reversals per Session');
    title('Average Reversals per Session by Age');
    grid on;
    
    disp('Plotted average reversals and reversals per session by age.');
end

