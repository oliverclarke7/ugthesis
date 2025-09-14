%% Reversal Experiment with Three Options (Renamed NP2, NP3, NP4)
% Session ends when 151 rewards are delivered.
% Up to exactly 4 reversals will occur.
% Each trial, the rat chooses one option; cumulative choices are recorded.
% A reward is delivered only if the chosen option is the current high reinforcement.
% A reversal is triggered when, over a rolling window of 30 trials, the rat chooses
% the current high option at least 21 times (when not in a latency period).
% At each reversal the high option is re-assigned and a latency period ensues.
% Each reversal is annotated on the plot with its type (e.g., "NP2 → NP3").

%% Parameters
maxRewards           = 151;    % End session after 151 rewards delivered
rollWindow           = 30;     % Rolling window size (in trials)
rollThreshold        = 21;     % Reversal triggered if high option is chosen >= 21 times in window
p_highChoice         = 0.65;   % Probability of choosing the current high option under normal conditions

% Post-reversal latency parameters
latencyTrials        = 25;     % Number of trials during the latency period
latencyNewHighProb   = 0.05;   % Initial probability of choosing the new high option
latencyOldHighProb   = 0.5;    % Initial probability of perseverating on old high option
adaptiveLatency      = true;   % Gradually adjust probabilities during latency

maxReversals         = 4;      % Allow up to 4 reversals

%% Initialization
currentHigh = 1;             % Start with NP2 (Option 1) as the high-reinforced option
previousHigh = NaN;          % Track old high option after reversal
latencyCounter = 0;
recentChoices = [];

% Cumulative choice counters
cumChoice1 = 0;    % NP2
cumChoice2 = 0;    % NP3
cumChoice3 = 0;    % NP4

% Session counters
totalRewards = 0;
trialCount   = 0;
numReversals = 0;

% For plotting
trialRecord       = [];
cumChoice1Record  = [];
cumChoice2Record  = [];
cumChoice3Record  = [];
reversalTrials    = [];
reversalLabels    = {};
xOffset = 5;

%% Simulation Loop
while totalRewards < maxRewards
    trialCount = trialCount + 1;

    %% 1. Determine the Rat's Choice
    if latencyCounter > 0
        trialsIntoLatency = latencyTrials - latencyCounter + 1;

        if adaptiveLatency
            probNewHigh = latencyNewHighProb + ...
                          (p_highChoice - latencyNewHighProb) * (trialsIntoLatency / latencyTrials);
            probOldHigh = latencyOldHighProb * (1 - trialsIntoLatency / latencyTrials);
        else
            probNewHigh = latencyNewHighProb;
            probOldHigh = latencyOldHighProb;
        end

        r = rand;
        if r < probOldHigh && ~isnan(previousHigh)
            choice = previousHigh;
        elseif r < (probOldHigh + probNewHigh)
            choice = currentHigh;
        else
            otherOptions = setdiff(1:3, [currentHigh, previousHigh]);
            choice = otherOptions(randi(numel(otherOptions)));
        end
    else
        if rand < p_highChoice
            choice = currentHigh;
        else
            otherOptions = setdiff(1:3, currentHigh);
            choice = otherOptions(randi(numel(otherOptions)));
        end
    end

    %% 2. Update Cumulative Choice Counts
    switch choice
        case 1
            cumChoice1 = cumChoice1 + 1;
        case 2
            cumChoice2 = cumChoice2 + 1;
        case 3
            cumChoice3 = cumChoice3 + 1;
    end

    %% 3. Reward Delivery
    if choice == currentHigh
        totalRewards = totalRewards + 1;
        outcome = 1;
    else
        outcome = 0;
    end

    %% 4. Reversal Logic
    if (latencyCounter == 0) && (numReversals < maxReversals)
        recentChoices = [recentChoices, outcome];
        if length(recentChoices) > rollWindow
            recentChoices = recentChoices(2:end);
        end

        if (length(recentChoices) == rollWindow) && (sum(recentChoices) >= rollThreshold)
            reversalTrials(end+1) = trialCount;
            numReversals = numReversals + 1;

            oldOptionLabel = sprintf('NP%d', currentHigh + 1);
            possibleNewHigh = setdiff(1:3, currentHigh);
            newHigh = possibleNewHigh(randi(numel(possibleNewHigh)));
            newOptionLabel = sprintf('NP%d', newHigh + 1);
            reversalLabels{end+1} = sprintf('%s → %s', oldOptionLabel, newOptionLabel);

            previousHigh = currentHigh;
            currentHigh = newHigh;
            latencyCounter = latencyTrials;
            recentChoices = [];
        end
    else
        if latencyCounter > 0
            latencyCounter = latencyCounter - 1;
        end
    end

    %% 5. Record Data for Plotting
    trialRecord(end+1)       = trialCount;
    cumChoice1Record(end+1)  = cumChoice1;
    cumChoice2Record(end+1)  = cumChoice2;
    cumChoice3Record(end+1)  = cumChoice3;
end

%% Plot the Cumulative Choice Curves
figure;
hold on;

plot(trialRecord, cumChoice1Record, '-', 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'NP2');
plot(trialRecord, cumChoice2Record, '-', 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'NP3');
plot(trialRecord, cumChoice3Record, '-', 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'NP4');

xlabel('Trial Number');
ylabel('Cumulative Choices');
title('Reversal Experiment: Cumulative Choices (NP2, NP3, NP4)');
legend('show', 'Location', 'NorthWest');
grid on;

if ~isempty(reversalTrials)
    yLimits = ylim;
    yPos = yLimits(2) * 0.80;

    hRev = xline(reversalTrials(1), 'r--', 'LineWidth', 2, 'DisplayName', 'Reversal');
    text(reversalTrials(1)+xOffset, yPos, reversalLabels{1}, 'Color', 'r', 'FontWeight', 'bold', ...
         'Rotation', 90, 'VerticalAlignment', 'middle');

    for i = 2:length(reversalTrials)
        xline(reversalTrials(i), 'r--', 'LineWidth', 2, 'HandleVisibility', 'off');
        text(reversalTrials(i)+xOffset, yPos, reversalLabels{i}, 'Color', 'r', 'FontWeight', 'bold', ...
             'Rotation', 90, 'VerticalAlignment', 'middle');
    end
end

hold off;