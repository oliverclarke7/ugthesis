%% compare_choice_probability_and_glme.m
% Aligns choice data to reversals, plots mean ± SEM by age, 
% then fits a mixed-effects logistic regression and prints ANOVA.

% ===== USER PARAMETERS =====
agesToPlot     = [30, 50, 70, 90, 120, 190];                % Ages to include
segmentPattern = 'REV1_REV';                              % Substring to filter segments
preTrials      = 40;                                      % Trials before reversal
postTrials     = 40;                                      % Trials after reversal
showSEM        = true;                                    % Toggle shaded SEM
saveFig        = true;                                    % Toggle saving figure
savePath       = '/Users/olly/Desktop/ASI/ASI/MAT Files/plots/';
markerDir      = '/Users/olly/Desktop/ASI/ASI/MAT Files/segment4.15';
segmentDir     = markerDir;                              % same folder for .txt files

% ===== ALIGN DATA BY AGE =====
xWindow         = -preTrials : postTrials;
allAlignedByAge = cell(1, numel(agesToPlot));

for a = 1:numel(agesToPlot)
    age = agesToPlot(a);
    markerFile = fullfile(markerDir, sprintf('REV_MARKERS_P%d.txt', age));
    if ~isfile(markerFile)
        warning('Skipping age %d: missing marker file.', age);
        continue;
    end
    mdata = readtable(markerFile, 'Delimiter', 'tab');
    mask  = contains(mdata.Segment_File, sprintf('_P%d_',age)) & ...
            contains(mdata.Segment_File, segmentPattern);
    files = mdata.Segment_File(mask);
    aligned = [];

    for i = 1:numel(files)
        fname = files{i};
        fpath = fullfile(segmentDir, fname);
        if ~isfile(fpath), warning('Missing file: %s', fname); continue; end
        data = readmatrix(fpath);
        if isempty(data), warning('Empty file: %s', fname); continue; end

        idx = find(strcmp(mdata.Segment_File, fname), 1);
        rev = mdata.Relative_Middle_Reversal(idx);
        T   = numel(data);

        row = nan(1, numel(xWindow));
        for j = 1:numel(xWindow)
            ti = rev + xWindow(j);
            if ti >= 1 && ti <= T
                row(j) = data(ti);
            end
        end
        aligned = [aligned; row]; %#ok<AGROW>
    end

    allAlignedByAge{a} = aligned;
end

% ===== PLOT MEAN ± SEM CURVES =====
% Define one RGB row per age (values 0–1):
colors = [
   0.50, 0.50, 0.50;  % P30 – grey
   0.85, 0.33, 0.10;  % P50 – brick red
   0.93, 0.69, 0.13;  % P70 – gold
   0.49, 0.18, 0.56;  % P90 – purple
   0.00, 0.50, 0.00;  % P120 – green
   0.30, 0.75, 0.93;  % P190 – sky blue
];
% (Alternatively, use a built-in map:
% colors = parula(numel(agesToPlot));
% )

figure('Color','w'); hold on; box on;
for a = 1:numel(agesToPlot)
    aligned = allAlignedByAge{a};
    if isempty(aligned), continue; end
    m = nanmean(aligned,1);
    s = nanstd(aligned,0,1) ./ sqrt(sum(~isnan(aligned),1));

    if showSEM
        % suppress shaded patch in legend:
        fill([xWindow, fliplr(xWindow)], ...
             [m+s, fliplr(m-s)], ...
             colors(a,:), ...
             'FaceAlpha', 0.2, ...
             'EdgeColor','none', ...
             'HandleVisibility','off');
    end

    plot(xWindow, m, ...
         'Color', colors(a,:), ...
         'LineWidth', 2, ...
         'DisplayName', sprintf('P%d', agesToPlot(a)));
end

% suppress reference line in legend
xline(0, 'k:', 'LineWidth', 1.5, 'HandleVisibility','off');

xlabel('Trial from reversal');
ylabel('p(NP_{best})');
title(sprintf('Choice of Highest NP aligned to %s', segmentPattern));
legend('Location','southeast');
ylim([0 1]);
xlim([min(xWindow) max(xWindow)]);
grid on;

if saveFig
    if ~exist(savePath, 'dir'), mkdir(savePath); end
    figName = fullfile(savePath, ...
        sprintf('choice_prob_aligned_to_reversal_%s.png', segmentPattern));
    saveas(gcf, figName);
    fprintf('Figure saved to: %s\n', figName);
end

% ===== BUILD LONG TABLE FOR GLME =====
Age     = [];
Session = [];
Time    = [];
Choice  = [];
sessID  = 0;

for a = 1:numel(agesToPlot)
    aligned = allAlignedByAge{a};
    [nSess, nTime] = size(aligned);
    for s = 1:nSess
        sessID = sessID + 1;
        for t = 1:nTime
            Age(end+1,1)     = agesToPlot(a);   %#ok<SAGROW>
            Session(end+1,1) = sessID;           %#ok<SAGROW>
            Time(end+1,1)    = xWindow(t);       %#ok<SAGROW>
            Choice(end+1,1)  = aligned(s,t);     %#ok<SAGROW>
        end
    end
end

tbl = table(Age, Session, Time, Choice);

% ===== FIT AND TEST MIXED-EFFECTS MODEL =====
glme = fitglme(tbl, ...
    'Choice ~ Age*Time + (1|Session)', ...
    'Distribution','Binomial', 'Link','logit');

disp('ANOVA for fixed effects (DFMethod = residual):');
disp(anova(glme,'DFMethod','residual'));