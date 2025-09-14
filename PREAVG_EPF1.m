%% explore_find_1_withLME_andPostHoc.m
% Extends explore_find_1.m by adding:
%  • Linear mixed-effects ANOVA for explore latency across ages
%  • Overall mean±SEM printing by age
%  • Shaded‐SEM ribbon plot in blue with close mean labels
%  • Paired post-hoc t-tests across all age pairs

%% PARAMETERS
critVal                   = 1.5;                  % critical value for file
maxTrialsBeforeReversal   = 10;                   % trials before reversal window
exploreSlopeThreshold     = [];                   % optional slope threshold
groupsToInclude           = {'SOC'};              % group filter (empty => all)

dataDir  = '/Users/olly/Desktop/ASI/ASI/MAT Files/';
dataFile = sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal);

%% 1) LOAD & FILTER
data = readtable(fullfile(dataDir, dataFile));
data.Properties.VariableNames = ...
    cellfun(@lower, data.Properties.VariableNames,'UniformOutput',false);
if ~isempty(groupsToInclude)
    data = data(ismember(data.group, groupsToInclude), :);
end

%% 2) COMPUTE EXPLORE LATENCY
nSess = height(data);
exploreLatency = nan(nSess,1);
for i = 1:nSess
    revTr = data.reversal_trial(i);

    % gather all CPs
    trials_before = [data.cp1_trial_before(i), data.cp2_trial_before(i), data.cp3_trial_before(i), data.cp4_trial_before(i)];
    slopes_before = [data.cp1_slope_before(i), data.cp2_slope_before(i), data.cp3_slope_before(i), data.cp4_slope_before(i)];
    trials_after  = [data.cp1_trial_after(i),  data.cp2_trial_after(i),  data.cp3_trial_after(i),  data.cp4_trial_after(i)];
    slopes_after  = [data.cp1_slope_after(i),  data.cp2_slope_after(i),  data.cp3_slope_after(i),  data.cp4_slope_after(i)];

    allTr = [trials_before, trials_after];
    allSl = [slopes_before,  slopes_after];
    valid = ~isnan(allTr) & ~isnan(allSl);
    allTr = allTr(valid);  allSl = allSl(valid);

    % sort by trial
    [allTr, si] = sort(allTr); allSl = allSl(si);

    % apply window
    wstart = revTr - maxTrialsBeforeReversal;
    idxW   = allTr >= wstart;
    allTr   = allTr(idxW);
    allSl   = allSl(idxW);

    % detect first drop
    cand = NaN;
    for j = 2:numel(allTr)
        if allSl(j) < allSl(j-1)
            if ~isempty(exploreSlopeThreshold) && ~isinf(exploreSlopeThreshold)
                if allSl(j) >= exploreSlopeThreshold
                    continue;
                end
            end
            cand = allTr(j) - revTr;
            break;
        end
    end
    exploreLatency(i) = cand;
end
data.explore_latency = exploreLatency;

%% 3) WITHIN-RAT AVERAGING (filtering out outliers)
[G, ratID, ratAge, ratGroup] = findgroups(data.rat, data.age, data.group);
ratExplore = splitapply(@(x) mean(rmoutliers(x), 'omitnan'), data.explore_latency, G);
T_rat = table(ratID, ratAge, ratGroup, ratExplore, ...
              'VariableNames', {'Rat','Age','Group','Latency'});

%% 4) PRINT OVERALL MEAN±SEM BY AGE
ages    = unique(T_rat.Age);
[Ga,~]  = findgroups(T_rat.Age);
meanLat = splitapply(@(x) mean(x,'omitnan'), T_rat.Latency, Ga);
semLat  = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T_rat.Latency, Ga);

fprintf('\nOverall Mean ± SEM by Age (Explore Latency):\n');
for k = 1:numel(ages)
    fprintf(' Age %3d: %.2f ± %.2f trials\n', ages(k), meanLat(k), semLat(k));
end

%% 5) PRETTIER SHADED-SEM RIBBON PLOT
x         = ages;
y         = meanLat;
e         = semLat;
lineColor = [0.2 0.6 0.8];
alphaVal  = 0.2;  % Less transparent now
dx        = 1;    % small horizontal label offset

figure('Color','w','Position',[200 200 600 400]);
hold on;

% ribbon (blue, now more visible)
fill([x; flipud(x)], [y-e; flipud(y+e)], ...
     lineColor, 'LineStyle','none', 'FaceAlpha', alphaVal);

% main line + markers
plot(x, y, '-o', ...
     'LineWidth',2, ...
     'MarkerSize',8, ...
     'MarkerFaceColor', lineColor, ...
     'Color', lineColor);

% tight mean labels (20% of SEM above point)
for i = 1:numel(x)
    dy = e(i)*0.2;
    text(x(i)+dx, y(i)+dy, sprintf('%.2f', y(i)), ...
         'FontSize',10, 'Color',[0.2 0.2 0.2], ...
         'HorizontalAlignment','left', 'VerticalAlignment','bottom');
end

% axes styling
ax = gca;
ax.Box       = 'off';
ax.TickDir   = 'out';
ax.LineWidth = 1.2;
ax.FontSize  = 12;
ax.FontName  = 'Helvetica';
ax.XTick     = x;
ax.YLim      = [0, max(y+e)*1.1];

xlabel('Age (days)',             'FontSize',14,'FontWeight','normal');
ylabel('Explore Latency (trials)','FontSize',14,'FontWeight','normal');
title('CP Explore Latency',      'FontSize',16,'FontWeight','normal');

grid on;
ax.GridColor = [0.8 0.8 0.8];
ax.GridAlpha = 0.5;

hold off;

%% 6) LME ANOVA ACROSS AGES
T_lme       = T_rat;
T_lme.Rat   = categorical(T_lme.Rat);
T_lme.Age_f = categorical(T_lme.Age);

fprintf('\n--- LME ANOVA: Explore Latency ~ Age ---\n');
lme = fitlme(T_lme, 'Latency ~ Age_f + (1|Rat)');
disp(anova(lme));

%% 7) PAIRED POST-HOC T-TESTS ACROSS AGE PAIRS
pairs = nchoosek(ages,2);
fprintf('\n--- Paired t-tests: Explore Latency ---\n');
for i = 1:size(pairs,1)
    a1 = pairs(i,1); a2 = pairs(i,2);
    r1 = T_rat.Rat(T_rat.Age==a1);
    r2 = T_rat.Rat(T_rat.Age==a2);
    C  = intersect(r1, r2);
    x  = T_rat.Latency(ismember(T_rat.Rat,C) & T_rat.Age==a1);
    y  = T_rat.Latency(ismember(T_rat.Rat,C) & T_rat.Age==a2);
    [~,p,~,st] = ttest(x,y);
    fprintf(' %3d vs %3d: t=%.2f, df=%d, p=%.3f\n', a1, a2, st.tstat, st.df, p);
end