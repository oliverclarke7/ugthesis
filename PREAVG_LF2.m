%% learn_find_2_withLME_andPostHoc.m
% Extends learn_find_2.m by adding:
%  • Printing mean ± SEM for each age and metric
%  • Linear mixed-effects ANOVA for each CP-before metric
%  • Paired post-hoc t-tests across ages for each metric
%  • Shaded-SEM “ribbon” plots using trapezoids per segment,
%    with ribbon color matching the line (but very transparent)

%% PARAMETERS
critVal                 = 1.5;            
slopeThreshold          = 0.6;            
minTrialsBeforeNextCP   = 5;              
minTrialsBeforeReversal = 0;              
maxTrialsBeforeReversal = 35;             
reversalFilter          = 'REV1_REV2';    
groupsToInclude         = {'SOC'};        

%% 1) LOAD & FILTER
dataDir  = '/Users/olly/Desktop/ASI/ASI/MAT Files/';
dataFile = sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal);
data     = readtable(fullfile(dataDir, dataFile));

% convert variable names to lowercase
data.Properties.VariableNames = ...
    cellfun(@lower, data.Properties.VariableNames, 'UniformOutput', false);

data     = data(strcmp(data.reversal_name, reversalFilter), :);
data     = data(ismember(data.group, groupsToInclude), :);

%% 2) SESSION-LEVEL CP-BEFORE METRICS
n = height(data);
rawFirstCandidateSlope     = nan(n,1);
weightedCandidateSlope     = nan(n,1);
weightedCandidateSlopeFig3 = nan(n,1);

for i = 1:n
    slopes = [data.cp1_slope_before(i), data.cp2_slope_before(i), ...
              data.cp3_slope_before(i), data.cp4_slope_before(i)];
    trials = [data.cp1_trial_before(i), data.cp2_trial_before(i), ...
              data.cp3_trial_before(i), data.cp4_trial_before(i)];
    revTr  = data.reversal_trial(i);
    
    candSl = []; candWt = []; firstFound = false;
    for j = 1:4
        s = slopes(j); t = trials(j);
        if isnan(s)|| isnan(t) || s < slopeThreshold, continue; end
        gap = revTr - t;
        if gap < minTrialsBeforeReversal || gap > maxTrialsBeforeReversal, continue; end
        if j<4 && ~isnan(trials(j+1)), nextEvt = trials(j+1); else nextEvt = revTr; end
        if nextEvt - t < minTrialsBeforeNextCP, continue; end
        
        candSl(end+1) = s; %#ok<AGROW>
        candWt(end+1) = nextEvt - t; %#ok<AGROW>
        if ~firstFound
            rawFirstCandidateSlope(i) = s;
            firstFound = true;
        end
    end
    if ~isempty(candSl)
        weightedCandidateSlope(i) = sum(candSl .* candWt) / sum(candWt);
    end
    
    if ~isnan(rawFirstCandidateSlope(i))
        fcIdx = find(slopes == rawFirstCandidateSlope(i), 1);
        allSl = []; allWt = [];
        for j = fcIdx:4
            t = trials(j);
            if isnan(t) || t >= revTr, continue; end
            s = slopes(j);
            if j<4 && ~isnan(trials(j+1)) && trials(j+1) < revTr
                w = trials(j+1) - t;
            else
                w = revTr - t;
            end
            allSl(end+1) = s; allWt(end+1) = w; %#ok<AGROW>
        end
        if ~isempty(allSl)
            weightedCandidateSlopeFig3(i) = sum(allSl .* allWt) / sum(allWt);
        end
    end
end

data.rawFirst  = rawFirstCandidateSlope;
data.weighted1 = weightedCandidateSlope;
data.weighted2 = weightedCandidateSlopeFig3;

%% 3) WITHIN-RAT AVERAGING
[G, ratID, ratAge, ratGroup] = findgroups(data.rat, data.age, data.group);
ratRaw       = splitapply(@(x) mean(x,'omitnan'), data.rawFirst, G);
ratWeighted1 = splitapply(@(x) mean(x,'omitnan'), data.weighted1, G);
ratWeighted2 = splitapply(@(x) mean(x,'omitnan'), data.weighted2, G);

T_rat = table(ratID, ratAge, ratGroup, ratRaw, ratWeighted1, ratWeighted2, ...
              'VariableNames', {'Rat','Age','Group','Raw','W1','W2'});

%% 4) AGE-LEVEL MEANS & SEM (for printing & plotting)
[GA, ages] = findgroups(T_rat.Age);
meanRaw = splitapply(@(x) mean(x,'omitnan'),            T_rat.Raw,  GA);
semRaw  = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T_rat.Raw,  GA);
meanW1  = splitapply(@(x) mean(x,'omitnan'),            T_rat.W1,   GA);
semW1   = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T_rat.W1,   GA);
meanW2  = splitapply(@(x) mean(x,'omitnan'),            T_rat.W2,   GA);
semW2   = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T_rat.W2,   GA);

%% 5) PRETTIER SHADED-SEM RIBBON (per-segment trapezoids,
%   ribbon color matches lineColor, FaceAlpha=0.05)
titles = { ...
    'Raw First Candidate Slope', ...
    'Weighted Candidate Slope', ...
    'Gross Weighted Slope (All CPs)' };
means = { meanRaw, meanW1, meanW2 };
sems  = { semRaw,  semW1,  semW2 };
lineColor = [0.2 0.6 0.8];

for m = 1:3
    x = ages;
    y = means{m};
    e = sems{m};
    
    figure('Color','w','Position',[200 200 600 400]);
    hold on;
    
    % trapezoid ribbon segments
    for k = 1:(numel(x)-1)
        xx = [x(k), x(k+1), x(k+1), x(k)];
        yy = [y(k)+e(k), y(k+1)+e(k+1), y(k+1)-e(k+1), y(k)-e(k)];
        fill(xx, yy, lineColor, 'EdgeColor','none', 'FaceAlpha',0.05);
    end
    
    % main line + markers
    plot(x, y, '-o', ...
         'LineWidth',2, ...
         'MarkerSize',8, ...
         'MarkerFaceColor',lineColor, ...
         'Color',lineColor);
    
    % value labels
    dx = 2; dy = max(e)*0.05;
    for i = 1:numel(x)
        if ~isnan(y(i))
            text(x(i)+dx, y(i)+e(i)+dy, sprintf('%.2f', y(i)), ...
                 'FontSize',10, 'Color',[0.2 0.2 0.2], ...
                 'HorizontalAlignment','left', 'VerticalAlignment','bottom');
        end
    end
    
    % axes styling
    ax = gca;
    ax.Box       = 'off';
    ax.TickDir   = 'out';
    ax.LineWidth = 1.2;
    ax.FontSize  = 12;
    ax.FontName  = 'Helvetica';
    ax.XTick     = x;
    ax.YLim      = [min(y-e)*0.9, max(y+e)*1.1];
    
    xlabel('Age (days)',    'FontSize',14,'FontWeight','normal');
    ylabel('Slope',         'FontSize',14,'FontWeight','normal');
    title(titles{m},       'FontSize',16,'FontWeight','normal');
    
    grid on;
    ax.GridColor = [0.8 0.8 0.8];
    ax.GridAlpha = 0.5;
    
    hold off;
end

%% 6) LINEAR MIXED-EFFECTS ANOVA
T_lme       = T_rat;
T_lme.Rat   = categorical(T_lme.Rat);
T_lme.Age_f = categorical(T_lme.Age);

fprintf('\n--- LME ANOVA: Raw ---\n');
disp(anova(fitlme(T_lme, 'Raw ~ Age_f + (1|Rat)')));
fprintf('\n--- LME ANOVA: W1 ---\n');
disp(anova(fitlme(T_lme, 'W1  ~ Age_f + (1|Rat)')));
fprintf('\n--- LME ANOVA: W2 ---\n');
disp(anova(fitlme(T_lme, 'W2  ~ Age_f + (1|Rat)')));

%% 7) POST-HOC PAIRED T-TESTS ACROSS AGES
pairs = nchoosek(ages,2);
runPaired(T_rat, pairs, 'Raw', 'Raw First Candidate Slope');
runPaired(T_rat, pairs, 'W1',  'Weighted Candidate Slope');
runPaired(T_rat, pairs, 'W2',  'Gross Weighted Slope');

%% Local helper: paired t-tests by rat
function runPaired(tbl, pairs, field, name)
    fprintf('\n--- Post-hoc paired t-tests: %s ---\n', name);
    for k = 1:size(pairs,1)
        a1 = pairs(k,1); a2 = pairs(k,2);
        common = intersect(tbl.Rat(tbl.Age==a1), tbl.Rat(tbl.Age==a2));
        x = tbl.(field)(ismember(tbl.Rat,common)&tbl.Age==a1);
        y = tbl.(field)(ismember(tbl.Rat,common)&tbl.Age==a2);
        [~,p,~,st] = ttest(x,y);
        fprintf('Age %d vs %d: t=%.2f, p=%.3f\n', a1, a2, st.tstat, p);
    end
end