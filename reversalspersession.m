%% plot_and_analyze_avg_reversals.m
% Loads per‐rat average reversals, filters by rat/group/sex, 
% prints mean & SEM, plots mean ± SEM across ages as a ribbon, 
% then runs a repeated‐measures ANOVA on Age and Bonferroni‐corrected paired t‐tests.

clearvars; close all; clc;

%% --- LOAD & CONFIGURATION ---

% Per‐rat reversal averages
data = readtable('/Users/olly/Desktop/ASI/ASI/MAT Files/avg_reversals_per_rat.csv');

% Rats to exclude
excluded_rats   = {'DEV019020016','DEV005017010'};

% Filters
included_groups = {'SOC'};          % {'ISO','SOC','SOC/ISO'}
included_sexes  = {'M','F'};        % {'M','F'}
included_ages   = [30 50 70 90 120 190];

% Core data for sex lookup
core       = readtable('/Users/olly/Desktop/ASI/ASI/MAT Files/core_ASI_data.csv');
sex_lookup = containers.Map(core.Rat, core.Sex);

%% --- FILTER DATA ---

% 1) Exclude unwanted rats
data( ismember(data.Rat, excluded_rats), : ) = [];

% 2) Keep only selected groups
data = data( ismember(data.Group, included_groups), : );

% 3) Keep only selected sexes
keep = false(height(data),1);
for i = 1:height(data)
    r = data.Rat{i};
    if isKey(sex_lookup, r) && ismember(sex_lookup(r), included_sexes)
        keep(i) = true;
    end
end
data = data(keep, :);

%% --- MERGE SEX INTO TABLE ---

sex_col = cell(height(data),1);
for i = 1:height(data)
    r = data.Rat{i};
    if isKey(sex_lookup, r)
        sex_col{i} = sex_lookup(r);
    else
        sex_col{i} = 'U';
    end
end
data.Sex = categorical(sex_col);

%% --- COMPUTE MEAN & SEM FOR PLOTTING ---

nAges    = numel(included_ages);
avg_vals = nan(1,nAges);
sem_vals = nan(1,nAges);

for k = 1:nAges
    age = included_ages(k);
    col = "AvgReversals_P" + string(age);
    if ismember(col, data.Properties.VariableNames)
        vals = data.(col);
        vals = vals(~isnan(vals));
        if ~isempty(vals)
            avg_vals(k) = mean(vals);
            sem_vals(k) = std(vals)/sqrt(numel(vals));
        end
    end
end

%% --- PRINT MEAN & SEM TABLE ---

fprintf('\nAge\tMean\tSEM\n');
for k = 1:nAges
    fprintf('%3d\t%0.2f\t%0.2f\n', included_ages(k), avg_vals(k), sem_vals(k));
end

%% --- PLOT MEAN ± SEM ACROSS AGES (with Ribbon) ---

figure('Color','w'); hold on; box on;

% Ribbon (shaded error region)
x = included_ages;
y = avg_vals;
y_upper = y + sem_vals;
y_lower = y - sem_vals;

fill([x fliplr(x)], [y_upper fliplr(y_lower)], [1 0 0], ... % Red color
    'FaceAlpha', 0.2, 'EdgeColor','none');                  % Transparent ribbon

% Plot mean line
plot(x, y, 'r-o', 'LineWidth', 2, 'MarkerSize', 6);

xlabel('Age (days)');
ylabel('Avg Reversals per Session');
title('Average Reversals Across Age (with SEM Ribbon)');
grid on;

% Annotate mean values
for k = 1:nAges
    if ~isnan(avg_vals(k))
        text(included_ages(k), avg_vals(k)-0.05, sprintf('%.2f',avg_vals(k)), ...
            'HorizontalAlignment','left','VerticalAlignment','top');
    end
end

%% --- REPEATED‐MEASURES ANOVA ON AGE ---

% Prepare wide‐format table
vars    = "AvgReversals_P" + string(included_ages);
tblWide = data(:, [{'Rat','Group','Sex'}, cellstr(vars)]);

% Define within‐subject factor
Within = table( categorical("P" + string(included_ages))', ...
                'VariableNames',{'Age'} );

% Fit repeated‐measures model
rm = fitrm(tblWide, sprintf('%s-%s ~ 1', vars(1), vars(end)), ...
           'WithinDesign', Within);

% Run ANOVA
ranovatbl = ranova(rm, 'WithinModel','Age');
disp('=== Repeated‐Measures ANOVA on Avg Reversals ===');
disp(ranovatbl);

%% --- POST‐HOC PAIRED T‐TESTS (BONFERRONI) ---

pairs   = nchoosek(1:nAges,2);
alpha_b = 0.05/size(pairs,1);
fprintf('\n=== Post‐hoc Paired t‐Tests (Bonferroni α = %.4f) ===\n', alpha_b);

for p = 1:size(pairs,1)
    i = pairs(p,1);
    j = pairs(p,2);
    c1 = vars(i); 
    c2 = vars(j);
    valid = ~isnan(data.(c1)) & ~isnan(data.(c2));
    diff  = data.(c1)(valid) - data.(c2)(valid);
    [h,pval,ci,stats] = ttest(diff, 0, 'Alpha', alpha_b);
    fprintf('  %s vs %s: t(%d)=%.2f, p=%.4f\n', ...
        c1, c2, stats.df, stats.tstat, pval);
end