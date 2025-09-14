%% exploit_find_1_withLME_andPostHoc.m
% Computes CP_before latency for EXPLOIT sessions, averages by rat & age,
% plots mean ± SEM with a matching-color shaded-SEM ribbon,
% runs a linear mixed-effects ANOVA across ages,
% and performs paired post-hoc t-tests across all age pairs.

%% PARAMETERS
critVal                 = 1.5;
slopeThreshold          = 0.6;
minTrialsBeforeNextCP   = 0;
minTrialsBeforeReversal = 0;
maxTrialsBeforeReversal = 35;   % limit candidate CPs to within 35 trials before reversal
excludedReversal        = 'REV1_REV2';
groupsToInclude         = {'SOC'};
dataDir                 = '/Users/olly/Desktop/ASI/ASI/MAT Files/';
dataFile                = sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal);

%% 1) LOAD & FILTER
data = readtable(fullfile(dataDir, dataFile));

% lowercase variable names
data.Properties.VariableNames = ...
    cellfun(@lower, data.Properties.VariableNames, 'UniformOutput', false);

% exclude specified reversal and non-SOC group
data = data(~strcmp(data.reversal_name, excludedReversal), :);
data = data(ismember(data.group, groupsToInclude), :);

%% 2) SESSION-LEVEL CP_before LATENCY
n = height(data);
cpBeforeLatency = nan(n,1);
for i = 1:n
    slopes  = [data.cp1_slope_before(i), data.cp2_slope_before(i), ...
               data.cp3_slope_before(i), data.cp4_slope_before(i)];
    trials  = [data.cp1_trial_before(i), data.cp2_trial_before(i), ...
               data.cp3_trial_before(i), data.cp4_trial_before(i)];
    revTr   = data.reversal_trial(i);
    candIdx = find(slopes >= slopeThreshold & ~isnan(slopes));
    for j = candIdx
        gapRev = revTr - trials(j);
        if gapRev < minTrialsBeforeReversal || gapRev > maxTrialsBeforeReversal
            continue
        end
        if j < numel(trials) && ~isnan(trials(j+1))
            nextEvt = trials(j+1);
        else
            nextEvt = revTr;
        end
        if (nextEvt - trials(j)) < minTrialsBeforeNextCP
            continue
        end
        cpBeforeLatency(i) = trials(j);
        break
    end
end
data.cp_before_latency = cpBeforeLatency;

%% 3) WITHIN-RAT AVERAGING (rat × age)
[G, ratID, ratAge, ratGroup] = findgroups(data.rat, data.age, data.group);

ratLatency = splitapply(@(x) mean(rmoutliers(x),'omitnan'), ...
                        data.cp_before_latency, G);

T_rat = table(ratID, ratAge, ratGroup, ratLatency, ...
              'VariableNames', {'RatID','Age','Group','Latency'});

%% 4) AGE-LEVEL MEAN & SEM
[GA, ages]    = findgroups(T_rat.Age);
meanLatency   = splitapply(@(x) mean(x,'omitnan'),      T_rat.Latency, GA);
semLatency    = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T_rat.Latency, GA);

fprintf('\nMean ± SEM by Age (CP_before latency):\n');
for k = 1:numel(ages)
    fprintf('  Age %d: Mean = %.2f, SEM = %.2f\n', ...
            ages(k), meanLatency(k), semLatency(k));
end

%% 5) PRETTIER SHADED-SEM RIBBON PLOT
x = ages;
y = meanLatency;
e = semLatency;
lineColor = [0.2 0.6 0.8];
alphaVal  = 0.2;  % More visible SEM ribbon

figure('Color','w','Position',[200 200 600 400]);
hold on;

% draw trapezoid ribbon segments with matching color and higher alpha
for k = 1:(numel(x)-1)
    xx = [x(k), x(k+1), x(k+1), x(k)];
    yy = [y(k)+e(k), y(k+1)+e(k+1), y(k+1)-e(k+1), y(k)-e(k)];
    fill(xx, yy, lineColor, 'EdgeColor','none', 'FaceAlpha', alphaVal);
end

% main line + markers
plot(x, y, '-o', ...
     'LineWidth',2, ...
     'MarkerSize',8, ...
     'MarkerFaceColor',lineColor, ...
     'Color',lineColor);

% value labels
dx = 2; dy = 1.0;
for i = 1:numel(x)
    if ~isnan(y(i))
        text(x(i)+dx, y(i)+dy, sprintf('%.2f',y(i)), ...
             'FontSize',10, 'Color',[0.2 0.2 0.2], ...
             'HorizontalAlignment','left', 'VerticalAlignment','bottom');
    end
end

% axis formatting
ax = gca;
ax.Box       = 'off';
ax.TickDir   = 'out';
ax.LineWidth = 1.2;
ax.FontSize  = 12;
ax.FontName  = 'Helvetica';
ax.XTick     = x;
ax.YLim      = [0, max(y+e)*1.1];

xlabel('Age (days)',       'FontSize',14,'FontWeight','normal');
ylabel('Latency (trials)', 'FontSize',14,'FontWeight','normal');
title('CP Exploit Latency – Excluding REV1', ...
      'FontSize',16,'FontWeight','normal');

grid on;
ax.GridColor = [0.8 0.8 0.8];
ax.GridAlpha = 0.5;

hold off;

%% 6) LINEAR MIXED-EFFECTS ANOVA ACROSS AGES
T_lme       = T_rat;
T_lme.RatID = categorical(T_lme.RatID);
T_lme.Age_f = categorical(T_lme.Age);
fprintf('\n--- LME ANOVA: CP_before latency by Age ---\n');
disp(anova(fitlme(T_lme, 'Latency ~ Age_f + (1|RatID)')));

%% 7) POST-HOC PAIRED T-TESTS ACROSS AGES
pairs = nchoosek(ages,2);
fprintf('\n--- Post-hoc paired t-tests (CP_before latency) ---\n');
for idx = 1:size(pairs,1)
    a1 = pairs(idx,1); a2 = pairs(idx,2);
    common = intersect(T_rat.RatID(T_rat.Age==a1), T_rat.RatID(T_rat.Age==a2));
    x = T_rat.Latency(ismember(T_rat.RatID,common)&T_rat.Age==a1);
    y = T_rat.Latency(ismember(T_rat.RatID,common)&T_rat.Age==a2);
    [~,p,~,st] = ttest(x,y);
    fprintf('  Age %d vs %d: t=%.2f, df=%d, p=%.3f\n', ...
            a1, a2, st.tstat, st.df, p);
end