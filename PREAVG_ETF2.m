%% exploit_find_2_withLME_andPostHoc.m
% Computes three CP-before metrics for EXPLOIT sessions, averages by rat & age,
% prints mean±SEM, plots each metric with blue shaded ribbons and tight mean labels,
% runs LME ANOVA across ages, and performs paired post-hoc t-tests across age pairs.

%% PARAMETERS
critVal                 = 1.5;            
slopeThreshold          = 0.6;            
minTrialsBeforeNextCP   = 5;              
minTrialsBeforeReversal = 0;              
maxTrialsBeforeReversal = 35;             
excludedReversal        = 'REV1_REV2';    
groupsToInclude         = {'SOC'};        

dataDir  = '/Users/olly/Desktop/ASI/ASI/MAT Files/';
dataFile = sprintf('exploit_explore_data_crit%.1f_filtered.csv', critVal);

%% 1) LOAD & FILTER
data = readtable(fullfile(dataDir, dataFile));
data.Properties.VariableNames = ...
    cellfun(@lower, data.Properties.VariableNames,'UniformOutput',false);
data = data(~strcmp(data.reversal_name, excludedReversal), :);
data = data(ismember(data.group, groupsToInclude), :);

%% 2) SESSION-LEVEL METRICS
nSess   = height(data);
rawSlope = nan(nSess,1);
wSlope1  = nan(nSess,1);
wSlope2  = nan(nSess,1);
for i = 1:nSess
    slopes = [data.cp1_slope_before(i), data.cp2_slope_before(i), ...
              data.cp3_slope_before(i), data.cp4_slope_before(i)];
    trials = [data.cp1_trial_before(i), data.cp2_trial_before(i), ...
              data.cp3_trial_before(i), data.cp4_trial_before(i)];
    revTr  = data.reversal_trial(i);

    % Raw and W1
    candS = []; candW = []; firstFound = false;
    for j = 1:4
        s = slopes(j); t = trials(j);
        if isnan(s)||isnan(t)|| s< slopeThreshold, continue; end
        gap = revTr - t;
        if gap<minTrialsBeforeReversal||gap>maxTrialsBeforeReversal, continue; end
        if j<4 && ~isnan(trials(j+1)), nextEvt = trials(j+1); else nextEvt = revTr; end
        if nextEvt - t < minTrialsBeforeNextCP, continue; end
        candS(end+1)=s; candW(end+1)=nextEvt-t; %#ok<AGROW>
        if ~firstFound, rawSlope(i)=s; firstFound=true; end
    end
    if ~isempty(candS)
        wSlope1(i)=sum(candS.*candW)/sum(candW);
    end

    % W2
    if ~isnan(rawSlope(i))
        fcIdx = find(slopes == rawSlope(i),1);
        allS=[]; allW=[];
        for j = fcIdx:4
            t = trials(j);
            if isnan(t)||t>=revTr, continue; end
            s = slopes(j);
            if j<4 && ~isnan(trials(j+1)) && trials(j+1)<revTr
                w = trials(j+1)-t; else w = revTr-t; end
            allS(end+1)=s; allW(end+1)=w; %#ok<AGROW>
        end
        if ~isempty(allS)
            wSlope2(i)=sum(allS.*allW)/sum(allW);
        end
    end
end

data.rawSlope = rawSlope;
data.wSlope1  = wSlope1;
data.wSlope2  = wSlope2;

%% 3) WITHIN-RAT AVERAGING
[G, ratID, ratAge, ratGrp] = findgroups(data.rat, data.age, data.group);
ratRaw = splitapply(@(x) mean(rmoutliers(x),'omitnan'), data.rawSlope, G);
ratW1  = splitapply(@(x) mean(rmoutliers(x),'omitnan'), data.wSlope1, G);
ratW2  = splitapply(@(x) mean(rmoutliers(x),'omitnan'), data.wSlope2, G);
T_rat  = table(ratID, ratAge, ratGrp, ratRaw, ratW1, ratW2, ...
               'VariableNames',{'Rat','Age','Group','Raw','W1','W2'});

%% 4) PRINT AGE-LEVEL MEANS & SEM
ages    = unique(T_rat.Age);
[G2,~]  = findgroups(T_rat.Age);
meanRaw = splitapply(@(x) mean(x,'omitnan'),      T_rat.Raw,  G2);
semRaw  = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T_rat.Raw,  G2);
meanW1  = splitapply(@(x) mean(x,'omitnan'),      T_rat.W1,   G2);
semW1   = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T_rat.W1,   G2);
meanW2  = splitapply(@(x) mean(x,'omitnan'),      T_rat.W2,   G2);
semW2   = splitapply(@(x) std(x,'omitnan')/sqrt(numel(x)), T_rat.W2,   G2);

fprintf('\nMean ± SEM by Age:\n');
for k = 1:numel(ages)
    fprintf(' Age %3d — Raw: %.2f±%.2f; W1: %.2f±%.2f; W2: %.2f±%.2f\n', ...
        ages(k), meanRaw(k), semRaw(k), meanW1(k), semW1(k), meanW2(k), semW2(k));
end

%% 5) PRETTIER PLOTS WITH BLUE RIBBONS & TIGHT MEAN LABELS
titles    = {'Raw First Candidate Slope', ...
             'Weighted Candidate Slope', ...
             'Gross Weighted Slope (All CPs)'};
means     = {meanRaw, meanW1, meanW2};
sems      = {semRaw, semW1, semW2};
lineColor = [0.2 0.6 0.8];
alphaVal  = 0.05;
dx        = 1;  % shift label slightly to the right

for m = 1:3
    x = ages; y = means{m}; e = sems{m};
    figure('Color','w','Position',[200 200 600 400]);
    hold on;
    % draw ribbon segments
    for k = 1:numel(x)-1
        xx = [x(k), x(k+1), x(k+1), x(k)];
        yy = [y(k)+e(k), y(k+1)+e(k+1), y(k+1)-e(k+1), y(k)-e(k)];
        fill(xx, yy, lineColor, 'EdgeColor','none', 'FaceAlpha',alphaVal);
    end
    % main line + markers
    plot(x, y, '-o', 'LineWidth',2, 'MarkerSize',8, ...
         'MarkerFaceColor',lineColor, 'Color',lineColor);
    % mean labels at 20% SEM above each point
    for i = 1:numel(x)
        dy = e(i)*0.2;
        text(x(i)+dx, y(i)+dy, sprintf('%.2f',y(i)), ...
             'FontSize',10, 'Color',[0.2 0.2 0.2], ...
             'HorizontalAlignment','left', 'VerticalAlignment','bottom');
    end
    xlabel('Age (days)','FontSize',14,'FontWeight','normal');
    ylabel('Slope','FontSize',14,'FontWeight','normal');
    title(titles{m},'FontSize',16,'FontWeight','normal');
    grid on;
    ax = gca;
    % ONLY for the Raw-slope plot, start y-axis at 0.82 and cap just above data
    if m == 1
        ax.YLim = [0.82, max(y+e)*1.02];
    end
    ax.Box       = 'off';
    ax.TickDir   = 'out';
    ax.LineWidth = 1.2;
    ax.FontSize  = 12;
    ax.FontName  = 'Helvetica';
    ax.XTick     = x;
    hold off;
end

%% 6) LME ANOVA ACROSS AGES
T_lme      = T_rat;
T_lme.Rat  = categorical(T_lme.Rat);
T_lme.AgeF = categorical(T_lme.Age);
fprintf('\n--- LME ANOVA: Raw ---\n');
disp(anova(fitlme(T_lme,'Raw ~ AgeF + (1|Rat)')));
fprintf('--- LME ANOVA: W1 ---\n');
disp(anova(fitlme(T_lme,'W1  ~ AgeF + (1|Rat)')));
fprintf('--- LME ANOVA: W2 ---\n');
disp(anova(fitlme(T_lme,'W2  ~ AgeF + (1|Rat)')));

%% 7) PAIRED POST-HOC T-TESTS ACROSS AGE PAIRS
pairs = nchoosek(ages,2);
for fld = {'Raw','W1','W2'}
    name = fld{1};
    fprintf('\nPaired t-tests: %s\n', name);
    for idx = 1:size(pairs,1)
        a1 = pairs(idx,1); a2 = pairs(idx,2);
        C  = intersect(T_rat.Rat(T_rat.Age==a1), T_rat.Rat(T_rat.Age==a2));
        x  = T_rat.(name)(ismember(T_rat.Rat,C)&T_rat.Age==a1);
        y  = T_rat.(name)(ismember(T_rat.Rat,C)&T_rat.Age==a2);
        [~,p,~,st] = ttest(x,y);
        fprintf(' %3d vs %3d: t=%.2f, df=%d, p=%.3f\n', a1,a2,st.tstat,st.df,p);
    end
end