%% pre_averages_learn_find_withLME_andPairedT.m
% Computes “CP before” latency, plots mean ± SEM by age with matching‐color ribbons,
% runs a linear mixed-effects ANOVA (Age as fixed, Rat as random),
% and performs paired t-tests comparing every pair of ages.

%% PARAMETERS
critVal                 = 1.5;            % critical value
slopeThreshold          = 0.6;            % slope cutoff for candidate CP
minTrialsBeforeNextCP   = 5;              % min gap to next CP
minTrialsBeforeReversal = 0;              % min gap to reversal
maxTrialsBeforeReversal = 35;             % max gap to reversal
reversalFilter          = 'REV1_REV2';    % only this reversal
groupsToInclude         = {'SOC'};        % only this group

%% 1) LOAD & FILTER
dataDir = '/Users/olly/Desktop/ASI/ASI/MAT Files/';
fname   = sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal);
T       = readtable(fullfile(dataDir, fname));

% lowercase variable names
T.Properties.VariableNames = ...
    cellfun(@lower, T.Properties.VariableNames,'UniformOutput',false);

T = T(strcmp(T.reversal_name, reversalFilter), :);
T = T(ismember(T.group, groupsToInclude), :);

%% 2) SESSION-LEVEL CP-BEFORE LATENCY
n = height(T);
cpLatency = nan(n,1);
for i = 1:n
    slopes = [T.cp1_slope_before(i), T.cp2_slope_before(i), ...
              T.cp3_slope_before(i), T.cp4_slope_before(i)];
    trials = [T.cp1_trial_before(i), T.cp2_trial_before(i), ...
              T.cp3_trial_before(i), T.cp4_trial_before(i)];
    revTr  = T.reversal_trial(i);
    cand   = find(slopes >= slopeThreshold & ~isnan(slopes));
    for k = cand
        gapRev = revTr - trials(k);
        if gapRev < minTrialsBeforeReversal || gapRev > maxTrialsBeforeReversal
            continue
        end
        if k < numel(trials) && ~isnan(trials(k+1))
            nextEvt = trials(k+1);
        else
            nextEvt = revTr;
        end
        if (nextEvt - trials(k)) < minTrialsBeforeNextCP
            continue
        end
        cpLatency(i) = trials(k);
        break
    end
end
T.cp_before_latency = cpLatency;

%% 3) RAT-LEVEL AVERAGES (per rat × age)
[G, ratID, ratAge] = findgroups(T.rat, T.age);
ratLatency = splitapply(@(x) mean(rmoutliers(x),'omitnan'), ...
                        T.cp_before_latency, G);
T_rat = table(ratID, ratAge, ratLatency, ...
              'VariableNames', {'Rat','Age','Latency'});

%% 4) AGE-LEVEL MEAN & SEM
[G2, ages] = findgroups(T_rat.Age);
meanLat = splitapply(@(x) mean(x,'omitnan'),      T_rat.Latency, G2);
semLat  = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T_rat.Latency, G2);

T_res = table(ages, meanLat, semLat, ...
   'VariableNames', {'Age','MeanLatency','SEM'});
disp(T_res);

%% 5) PRETTIER PLOT WITH MATCHING‐COLOR RIBBON
figure('Color','w','Position',[200 200 600 400]);
hold on;

% Define line/ribbon color
lineColor = [0.2 0.6 0.8];
alphaVal  = 0.2;  % More visible SEM ribbon
x = ages;
y = meanLat;
e = semLat;

% 1) shaded SEM ribbon
fill([x; flipud(x)], [y-e; flipud(y+e)], ...
     lineColor, 'LineStyle','none', 'FaceAlpha', alphaVal);

% 2) main line + markers
plot(x, y, '-o', ...
     'LineWidth',2, ...
     'MarkerSize',8, ...
     'MarkerFaceColor',lineColor, ...
     'Color',lineColor);

% 3) value labels just above each point
dx = 2;
dy = 1.0;  % shift labels up by ~1 trial
for i = 1:numel(x)
    if ~isnan(y(i))
        text(x(i)+dx, y(i)+dy, sprintf('%.2f', y(i)), ...
             'FontSize',10, 'Color',[0.2 0.2 0.2], ...
             'HorizontalAlignment','left', 'VerticalAlignment','bottom');
    end
end

% 4) axes formatting
ax = gca;
ax.Box       = 'off';
ax.TickDir   = 'out';
ax.LineWidth = 1.2;
ax.FontSize  = 12;
ax.FontName  = 'Helvetica';
ax.XTick     = x;
ax.YLim      = [0, max(y+e)*1.1];

% 5) labels & title
xlabel('Age (days)',       'FontSize',14, 'FontWeight','normal');
ylabel('Latency (trials)', 'FontSize',14, 'FontWeight','normal');
title('CP Learn Latency – REV1 Only', 'FontSize',16, 'FontWeight','normal');

% 6) light grid
grid on;
ax.GridColor = [0.8 0.8 0.8];
ax.GridAlpha = 0.5;

hold off;

%% 6) LINEAR MIXED-EFFECTS ANOVA
T_lme       = T_rat;
T_lme.Rat   = categorical(T_lme.Rat);
T_lme.Age_f = categorical(T_lme.Age);
lme         = fitlme(T_lme, 'Latency ~ Age_f + (1|Rat)');
anovaTbl    = anova(lme);
disp(anovaTbl);

%% 7) PAIRED T-TESTS BETWEEN AGES
T_rat.AgeStr = strcat('Age', string(T_rat.Age));
wideRat      = unstack(T_rat, 'Latency', 'AgeStr', 'GroupingVariable', 'Rat');
ageCols      = wideRat.Properties.VariableNames(2:end);
numAges      = numel(ageCols);

comparison = strings(nchoosek(numAges,2),1);
tstat      = zeros(size(comparison));
df         = zeros(size(comparison));
pval       = zeros(size(comparison));

idx = 1;
for i = 1:numAges-1
    for j = i+1:numAges
        x = wideRat.(ageCols{i});
        y = wideRat.(ageCols{j});
        [~,p,~,stats] = ttest(x,y);
        comparison(idx) = ageCols{i} + " vs " + ageCols{j};
        tstat(idx)      = stats.tstat;
        df(idx)         = stats.df;
        pval(idx)       = p;
        idx = idx + 1;
    end
end

pairedTbl = table(comparison, tstat, df, pval, ...
                   'VariableNames', {'Comparison','tStat','DF','pValue'});
disp(pairedTbl);