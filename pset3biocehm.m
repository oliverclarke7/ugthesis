% =========================================================================
% SCRIPT: compare_myoglobin_mutant.m
%
% DESCRIPTION:
%   This script plots the fractional saturation (Y) vs. partial pressure of 
%   oxygen (pO2) for both wild-type myoglobin and a distal His->Gly mutant.
%
%   The distal His->Gly mutation is presumed to lower the oxygen affinity 
%   (i.e., increase P50).
%
% =========================================================================

clc;        % Clear command window
clear;      % Clear variables

% -------------------------------------------------------------------------
% 1) Define range of oxygen partial pressures to display
%    (Units could be Torr, mmHg, or any consistent pressure unit)
% -------------------------------------------------------------------------
pO2 = linspace(0, 40, 500);  % from 0 to 40 in 500 increments

% -------------------------------------------------------------------------
% 2) Parameters for Wild-Type Myoglobin
%    - Typical P50 for myoglobin is ~2-3 Torr, but it can vary by source.
% -------------------------------------------------------------------------
P50_wt = 3.0;  % an example value (Torr)

% -------------------------------------------------------------------------
% 3) Parameters for Mutant (His->Gly)
%    - Loss of hydrogen bonding is expected to weaken O2 binding.
%    - This means the mutant has LOWER affinity, or equivalently HIGHER P50.
% -------------------------------------------------------------------------
P50_mut = 10.0;  % choose a higher value than wild type to simulate lower affinity

% -------------------------------------------------------------------------
% 4) Calculate fractional saturation (Y) for each
% -------------------------------------------------------------------------
Y_wt   = pO2 ./ (pO2 + P50_wt);
Y_mut  = pO2 ./ (pO2 + P50_mut);

% -------------------------------------------------------------------------
% 5) Plot the results
% -------------------------------------------------------------------------
figure('Color', 'w');  % new figure with white background

% Wild-type in blue
plot(pO2, Y_wt, 'b-', 'LineWidth', 2);
hold on;

% Mutant in red dashed
plot(pO2, Y_mut, 'r--', 'LineWidth', 2);

% -------------------------------------------------------------------------
% 6) Annotate the figure
% -------------------------------------------------------------------------
xlabel('Partial Pressure of O_2 (Torr)');
ylabel('Fractional Saturation, Y');
title('Myoglobin O_2 Binding: Wild-Type vs. His->Gly Mutant');
legend({'Wild-Type', 'Mutant'}, 'Location', 'southeast');
grid on;

% =========================================================================
% END OF SCRIPT
% =========================================================================