%% Analysis: Dataset after preprocessing
%
% version   1.0, 19.04.2024
% author    Miriam Kirchhoff
% project   C2B


%% Load data

path_current.load = '\simulation_data';
path_current.save = '\simulation_data';

load([path_current.load '\Summary_data_allsubjects'])

% Rejection stats total
fprintf(['Rejection: \n' ...
    '%.0f total \n' ...
    '%.0f EEG \n' ...
    '%.0f preinnervation \n' ...
    '%.0f no peak \n \n'], ...
    sum(all_results.rejection.total(results_particiant_idx)), ...
    sum(all_results.rejection.EEG_threshold(results_particiant_idx)), ...
    sum(all_results.rejection.preinnervation(results_particiant_idx)), ...
    sum(all_results.rejection.no_peak(results_particiant_idx)))

fprintf(['Rejection relative: \n' ...
    '%.2f %% total \n' ...
    '%.2f %% EEG \n' ...
    '%.2f %% preinnervation \n' ...
    '%.2f %% no peak \n \n'], ...
    sum(all_results.rejection.total(results_particiant_idx)         / sum(trial_particiant_idx) * 100), ...
    sum(all_results.rejection.EEG_threshold(results_particiant_idx) / sum(trial_particiant_idx) * 100), ...
    sum(all_results.rejection.preinnervation(results_particiant_idx)/ sum(trial_particiant_idx) * 100), ...
    sum(all_results.rejection.no_peak(results_particiant_idx)       / sum(trial_particiant_idx) * 100))


%% Analysis: Circular to linear regression

clear mdl_linear 
alpha = 0.05;

mep = zscore(all_trials.MEP_log(~all_trials.outlier_all & trial_particiant_idx));
phase = all_trials.phase(~all_trials.outlier_all & trial_particiant_idx);

% plot simulation data
mg_step = 0.05;                 % measurement grid stepsize
min_phase = -pi;                % min possible phase
max_phase = pi;                 % max possible phase
subject_data.transformed_mep_amplitudes = mep';
subject_data.phases                     = phase;
simulation_parameters.percentage_data = 0.2;
simulation_parameters.steps =   100;
simulation_parameters.dataset = [subject_data.phases',subject_data.transformed_mep_amplitudes];
simulation_parameters.wintype = 'winsize';
simulation_parameters.winsize = (max_phase-min_phase)*simulation_parameters.percentage_data;  % window size of simulation
simulation_parameters.limits =  [min_phase, max_phase];
shuffle_repetitions =           5;
[ground_truth, shuffle] = simulation_ground_truth_slidewindow(simulation_parameters, subject_data, shuffle_repetitions, 0, false, false);

mdl_linear.estphase = phase;    % enter phase values
mdl_linear.excitability = mep;     % enter MEPs
mdl_linear.fit = fitlm([cos(mdl_linear.estphase)' sin(mdl_linear.estphase)'], mdl_linear.excitability);

mdl_linear.pval = mdl_linear.fit.coefTest;
mdl_linear.R2ordinary = mdl_linear.fit.Rsquared.Ordinary;
mdl_linear.R2adjusted = mdl_linear.fit.Rsquared.Adjusted;

% plot regression line
% y = intercept + x1*cos(x) + x2 * sin(x)
x_reg = linspace(-pi, pi, 1000);
y_reg = mdl_linear.fit.Coefficients.Estimate(1) ...
    + mdl_linear.fit.Coefficients.Estimate(2).* cos(x_reg) ...
    + mdl_linear.fit.Coefficients.Estimate(3).* sin(x_reg);
mdl_linear.opt = x_reg(y_reg == max(y_reg));

% Disp results
disp('Group level model fit:')
if mdl_linear.fit.ModelFitVsNullModel.Pvalue < 0.001
    fprintf('R² = %.2f, F(%.0f, %.0f) = %.2f, p < 0.001 \n', ...
        mdl_linear.R2adjusted, ...
        mdl_linear.fit.NumPredictors, ...
        mdl_linear.fit.DFE, ...
        mdl_linear.fit.ModelFitVsNullModel.Fstat)
else
    fprintf('R² = %.2f, F(%.0f, %.0f) = %.2f, p = %.2f \n', ...
        mdl_linear.R2adjusted, ...
        mdl_linear.fit.NumPredictors, ...
        mdl_linear.fit.DFE, ...
        mdl_linear.fit.ModelFitVsNullModel.Fstat, ...
        mdl_linear.fit.ModelFitVsNullModel.Pvalue)
end
fprintf('Group level: Optimal phase = %.0f \n', rad2deg(mdl_linear.opt))


%% Single-subject analysis results

% get all phases and meps
phases = all_trials.phase(~all_trials.outlier_all & trial_particiant_idx);
meps = [];
for i = unique(all_trials.participant)
    meps = [meps, zscore(all_trials.MEP_log(all_trials.participant == i & ~all_trials.outlier_all & trial_particiant_idx))];
end
%meps = meps((all_trials.outlier_all<50));


%% Plotting

close all

clear colors
colors.lightblue    = [45 160 220]/255;
colors.darkgrey     = [1 1 1]*0.4;
colors.lightred     = [220 70 80]/255;
colors.lightgreen   = [220 220 40]/255;
colors.darkred      = [150 30 40]/255;
colors.darkblue     = [40 110 150]/255;
colors.lightgrey    = [1 1 1]*0.5;

subjects = unique(all_trials.participant(trial_particiant_idx));

% Benjamini-Hochberg
% is_significant = all_results.c2l_reg.p_corrected < alpha;

% Bonferroni
is_significant = all_results.c2l_reg.p_uncorrected*length(subjects) < alpha;

% Give significant subject results
disp('optimal phase')
disp(sort(rad2deg(all_results.c2l_reg.opt(is_significant & results_particiant_idx))))
disp('R²')
disp(sort(all_results.c2l_reg.R2ordinary(is_significant & results_particiant_idx)))
disp('mean and std R²')
disp(mean(all_results.c2l_reg.opt(is_significant & results_particiant_idx)))
disp(std(all_results.c2l_reg.opt(is_significant & results_particiant_idx)))



% Define plotting font
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'defaultLineLineWidth',2)


%% Paper figure: Group phase effect

ticks = linspace(-pi, pi, 9);

group_effect = figure('units','normalized','outerposition',[0.5 0 0.3 0.4]);
t = tiledlayout(1,2, 'TileSpacing','compact');
title(t, 'Relationship of MEP size and mu-phase', 'FontName', 'Times New Roman')
nexttile
grid on; hold on

% phase dependency across subjects
plot(ground_truth.X, ground_truth.mean, 'k', 'LineWidth', 2)
X_area = [ground_truth.X, flip(ground_truth.X)];
Y = [ground_truth.ci(1,:), flip(ground_truth.ci(2,:))];
fill(X_area, Y, 'm', 'LineStyle', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2)

% plot regression fit
plot(x_reg,y_reg, 'Color', colors.lightblue, 'LineWidth', 2)

% plot oscillation
x_osc = linspace(-pi, pi, 200);
y_osc = cos(x_osc)*range(Y)*0.1 + min(Y) - range(Y)*0.1 + 0.0028;
plot(x_osc, y_osc, 'color', [0 0 0] + 0.5, 'LineStyle', ':')
ylabel('MEP response (log, z-score) [µV]')
[rho, pval] = circ_corrcl(phases, meps);
legend(["Mean" "95% CI" "Regression fit" "Oscillation"], 'Location', 'southoutside', 'NumColumns',2)
xticks(ticks(1:end))
xticklabels({'trough','early rising','rising','late rising', 'peak','early falling','falling','late falling', 'trough'})
xtickangle(90)
xlim([-pi pi])
title('Group level', 'FontName', 'Times New Roman')
xlim([-pi pi])
ylim([-0.1 0.08])



% all subjects
nexttile; grid on; hold on
% optimal phase plot histogram
edges = ticks + 0.5*(ticks(2)- ticks(1));
edges = [-pi edges];
% plot points of optimal phase
scatter(200, 1, 'w', 'filled')
N_sub = size(all_results.c2l_reg.opt(results_particiant_idx));
scatter(all_results.c2l_reg.opt(results_particiant_idx)', (2*rand(N_sub)).*ones(N_sub), 'k')
x = linspace(-pi, pi, 200);
y = (cos(x) + 1)*0.5;
N = histcounts(all_results.c2l_reg.opt(results_particiant_idx),edges);

N([1 end]) = N(1) + N(end);
histogram('binedges', edges, 'bincounts', N, 'FaceColor','k', 'FaceAlpha', 0.1);
yticks([0:2:max(N)])
ylim([0 max(N)])
% plot(x, y, '.', 'color', [0 0 0] + 0.5)

% significant subjects
% optimal phase plot histogram
% plot points of optimal phase
scatter(200, 1, 'w', 'filled')
N_sub = size(all_results.c2l_reg.opt(is_significant & (results_particiant_idx)));
scatter(sort(all_results.c2l_reg.opt(is_significant & (results_particiant_idx)))', (2*rand(N_sub)).*ones(N_sub), 'k', 'filled')
x = linspace(-pi, pi, 200);
y = (cos(x) + 1)*0.5;
% histogram(all_results(is_significant,18),edges, 'FaceColor', 'k', 'FaceAlpha', 0.2)
N = histcounts(all_results.c2l_reg.opt(is_significant & (results_particiant_idx)),edges);

N([1 end]) = N(1) + N(end);
histogram('binedges', edges, 'bincounts', N, 'FaceColor',[1 1 1]*0.3);

xticks(ticks(1:end))
xticklabels({'trough','early rising','rising','late rising', 'peak','early falling','falling','late falling', 'trough'})
xtickangle(90)
ylabel('Number of subjects')
xlim([-pi pi])
legend(["\color{black}\bf All subjects" "Optimal phase" "Count" "\color{black}\bf Significant effect" "Optimal phase" "Count" "" ], 'location', 'southoutside', 'NumColumns',2)
title('Optimal phase location, all subjects', 'FontName', 'Times New Roman')

fontsize(gcf,scale=1.4)

group_effect.Color = 'w';

saveas(group_effect,[path_current.save '\figures\2phase.jpg'])

