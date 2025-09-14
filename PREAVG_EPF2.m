%% explore_find_2_withLME_andPostHoc.m
% Computes three explore CP slope metrics, averages by rat & age,
% prints mean±SEM, plots each metric with blue shaded ribbons and tight mean labels,
% runs LME ANOVA across ages, and performs paired post-hoc t-tests across age pairs.

%% PARAMETERS
critVal                  = 1.5;            % critical value for file selection
maxTrialsBeforeReversal  = 10;             % trials before reversal window
exploreSlopeThreshold    = [];             % optional slope threshold
groupsToInclude          = {'SOC'};        % group filter (empty ⇒ all)

dataDir  = '/Users/olly/Desktop/ASI/ASI/MAT Files/';
dataFile = sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal);

%% 1) LOAD & FILTER
data = readtable(fullfile(dataDir, dataFile));
data.Properties.VariableNames = ...
    cellfun(@lower, data.Properties.VariableNames,'UniformOutput',false);
if ~isempty(groupsToInclude)
    data = data(ismember(data.group,groupsToInclude), :);
end

ages      = unique(data.age);
nSessions = height(data);

%% 2) SESSION-LEVEL METRICS
absDecrease    = nan(nSessions,1);
pctDecrease    = nan(nSessions,1);
candidateSlope = nan(nSessions,1);

for i = 1:nSessions
    revTr = data.reversal_trial(i);
    tb    = [data.cp1_trial_before(i), data.cp2_trial_before(i), data.cp3_trial_before(i), data.cp4_trial_before(i)];
    sb    = [data.cp1_slope_before(i),  data.cp2_slope_before(i),  data.cp3_slope_before(i),  data.cp4_slope_before(i)];
    ta    = [data.cp1_trial_after(i),  data.cp2_trial_after(i),  data.cp3_trial_after(i),  data.cp4_trial_after(i)];
    sa    = [data.cp1_slope_after(i),   data.cp2_slope_after(i),   data.cp3_slope_after(i),   data.cp4_slope_after(i)];

    allTr = [tb, ta];
    allSl = [sb, sa];
    valid = ~isnan(allTr) & ~isnan(allSl);
    allTr = allTr(valid);
    allSl = allSl(valid);
    [allTr, si] = sort(allTr);
    allSl = allSl(si);

    % window
    ws   = revTr - maxTrialsBeforeReversal;
    mask = allTr >= ws;
    allTr = allTr(mask);
    allSl = allSl(mask);

    % detect first drop (negative change)
    a = NaN; p = NaN; s = NaN;
    for j = 2:numel(allTr)
        if allSl(j) < allSl(j-1)
            if ~isempty(exploreSlopeThreshold) && allSl(j) >= exploreSlopeThreshold
                continue;
            end
            a = allSl(j) - allSl(j-1);
            if allSl(j-1) > 0
                p = 100 * a / allSl(j-1);
            end
            s = allSl(j);
            break;
        end
    end

    absDecrease(i)    = a;
    pctDecrease(i)    = p;
    candidateSlope(i) = s;
end

data.abs_decrease    = absDecrease;
data.pct_decrease    = pctDecrease;
data.candidate_slope = candidateSlope;

%% 3) WITHIN-RAT AVERAGING
[G, ratID, ratAge, ratGroup] = findgroups(data.rat, data.age, data.group);
ratAbs  = splitapply(@(x) mean(rmoutliers(x),'omitnan'), data.abs_decrease,    G);
ratPct  = splitapply(@(x) mean(rmoutliers(x),'omitnan'), data.pct_decrease,    G);
ratCand = splitapply(@(x) mean(rmoutliers(x),'omitnan'), data.candidate_slope, G);

T_rat = table(ratID, ratAge, ratGroup, ratAbs, ratPct, ratCand, ...
              'VariableNames',{'Rat','Age','Group','AbsDec','PctDec','CandSlope'});

%% 4) PRINT MEAN ± SEM BY AGE
[Ga, ages] = findgroups(T_rat.Age);
meanAbs    = splitapply(@(x) mean(x,'omitnan'),               T_rat.AbsDec,   Ga);
semAbs     = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T_rat.AbsDec,   Ga);
meanPct    = splitapply(@(x) mean(x,'omitnan'),               T_rat.PctDec,   Ga);
semPct     = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T_rat.PctDec,   Ga);
meanCand   = splitapply(@(x) mean(x,'omitnan'),               T_rat.CandSlope,Ga);
semCand    = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T_rat.CandSlope,Ga);

fprintf('\nMean ± SEM by Age (Explore Metrics):\n');
for k = 1:numel(ages)
    fprintf(' Age %3d — AbsDec: %+0.2f±%0.2f; PctDec: %+0.2f±%0.2f; CandSlp: %.2f±%.2f\n', ...
        ages(k), meanAbs(k), semAbs(k), meanPct(k), semPct(k), meanCand(k), semCand(k));
end

%% 5) PRETTIER PLOTS WITH BLUE RIBBONS & TIGHT MEAN LABELS
metrics   = {meanAbs, meanPct, meanCand};
sems      = {semAbs, semPct, semCand};
titles    = {'Absolute Decrease in Slope', 'Percentage Decrease (%)', 'Candidate Slope After Decrease'};
lineColor = [0.2 0.6 0.8];
alphaVal  = 0.2;  % More visible SEM shading
dx        = 1;

for m = 1:3
    x = ages; y = metrics{m}; e = sems{m};

    figure('Color','w','Position',[200 200 600 400]); hold on;
    fill([x; flipud(x)], [y-e; flipud(y+e)], ...
         lineColor, 'EdgeColor','none', 'FaceAlpha', alphaVal);
    plot(x, y, '-o', 'LineWidth',2, 'MarkerSize',8, ...
         'MarkerFaceColor', lineColor, 'Color', lineColor);

    for i = 1:numel(x)
        dy = e(i)*0.2;
        if m < 3
            lbl = sprintf('%+0.2f', y(i));
        else
            lbl = sprintf('%.2f', y(i));
        end
        text(x(i)+dx, y(i)+dy, lbl, ...
             'FontSize',10, 'Color',[0.2 0.2 0.2], ...
             'HorizontalAlignment','left', 'VerticalAlignment','bottom');
    end

    ax = gca;
    ax.Box       = 'off';
    ax.TickDir   = 'out';
    ax.LineWidth = 1.2;
    ax.FontSize  = 12;
    ax.FontName  = 'Helvetica';
    ax.XTick     = x;

    switch m
        case 1  % absolute decrease
            ax.YLim = [min(y-e)*1.1, -0.35];
        case 2  % percentage decrease
            ax.YLim = [min(y-e)*1.1, -60];
        otherwise  % candidate slope
            ax.YLim = [min(y-e)*0.9, max(y+e)*1.1];
    end

    xlabel('Age (days)',             'FontSize',14,'FontWeight','normal');
    ylabel(titles{m},                'FontSize',14,'FontWeight','normal');
    title(titles{m},                 'FontSize',16,'FontWeight','normal');
    grid on; hold off;
end

%% 6) LME ANOVA FOR EACH METRIC
T_lme       = T_rat;
T_lme.Rat   = categorical(T_lme.Rat);
T_lme.Age_f = categorical(T_lme.Age);

fprintf('\n--- LME ANOVA: Absolute Decrease ---\n');
disp(anova(fitlme(T_lme,'AbsDec ~ Age_f + (1|Rat)')));
fprintf('--- LME ANOVA: Percentage Decrease ---\n');
disp(anova(fitlme(T_lme,'PctDec ~ Age_f + (1|Rat)')));
fprintf('--- LME ANOVA: Candidate Slope ---\n');
disp(anova(fitlme(T_lme,'CandSlope ~ Age_f + (1|Rat)')));

%% 7) PAIRED POST-HOC T-TESTS ACROSS AGE PAIRS
pairs = nchoosek(ages,2);
for fld = {'AbsDec','PctDec','CandSlope'}
    name = fld{1};
    fprintf('\nPaired t-tests: %s\n', name);
    for i = 1:size(pairs,1)
        a1 = pairs(i,1); a2 = pairs(i,2);
        C  = intersect(T_rat.Rat(T_rat.Age==a1), T_rat.Rat(T_rat.Age==a2));
        x  = T_rat.(name)(ismember(T_rat.Rat,C) & T_rat.Age==a1);
        y  = T_rat.(name)(ismember(T_rat.Rat,C) & T_rat.Age==a2);
        [~,p,~,st] = ttest(x,y);
        fprintf(' %3d vs %3d: t=%.2f, df=%d, p=%.3f\n', a1, a2, st.tstat, st.df, p);
    end
end