function [all_results, all_trials] = generate_example_data(settings)
% function [all_results, all_trials] = generate_example_data(settings)
%
% Generate example data for the simulation
%
% INPUTS
% settings:     Struct with fields
    % subjects: vector with subject indices (1 x n_subjects)
%
% OUTPUTS
% all_trials: struct with fields
    % phase         Vector, Phase of each trial (1 x n_trials)
    % MEP_log       Vector, MEP of each trial (1 x n_trials)
    % participant   Vector, participant nr of each trial (1 x n_trials)
    % outlier_all   Binary vector, outlier value of each trial (1 x n_trials)
% all_results: struct with fields
    % participant   Vector, participant nr (1 x n_participants)
    % moving_win.opt    Vector, optimal phase for moving window analysis (1 x n_participants)
    % c2l_reg.opt   Vector, optimal phase for regression (1 x n_participants)
    % c2l_reg.R2ordinary    Vector, R² for regression (1 x n_participants)
    % c2l_reg.p_uncorrected Vector, p-value for regression (1 x n_participants)
    % c2l_reg.p_corrected   Vector, corrected p-value for regression (1 x n_participants)
%
% version   1.0, 22.04.2024
% author    Miriam Kirchhoff
% project   C2B

% generate at least one subject
if isempty(settings.subjects)
    settings.subjects = 1;
end

% initialize output array
all_trials.phase = [];
all_trials.MEP_log = [];
all_trials.participant = [];
all_trials.outlier_all = [];

% number of samples per participant
N = 1200;

% loop over subjects
for idx_participant = 1:length(settings.subjects)

    mdl_linear.R2ordinary = 0;

    while mdl_linear.R2ordinary < 0.022 % effect size used in the paper
    % determine random parameters: phase shift, amplitude
    % Phase between +/- pi
    phase_shift = rand(1)*2*pi - pi;

    % Amplitude between 0.5 and 1
    amplitude = rand(1) * 0.5 + 0.5;

    % generate random phase locations
    phase = rand(1, N)*2*pi - pi;

    % generate MEP sizes
    mep_size = amplitude .* sin(phase - phase_shift);

    % add structural noise
    noise_phase_shift = rand(1)*2*pi - pi;
    noise_amplitude = rand(1) * 0.01 + 0.01;
    noise_struct = noise_amplitude .* sin(3*(phase) - noise_phase_shift);

    % add gaussian noise with mean 0 and sd 1 since data was z-scored
    noise_rand = normrnd(0,1, 1, N);
    mep_size = mep_size + noise_rand + noise_struct;
    participant = ones(size(phase)) * settings.subjects(idx_participant);

    % Analyse relationship between phase and MEP in circular to linear
    % regression
    mdl_linear.estphase = phase;    % enter phase values
    mdl_linear.excitability = mep_size;     % enter MEPs
    mdl_linear.fit = fitlm([cos(mdl_linear.estphase)' sin(mdl_linear.estphase)'], mdl_linear.excitability);
    mdl_linear.pval = mdl_linear.fit.coefTest;
    mdl_linear.R2ordinary = mdl_linear.fit.Rsquared.Ordinary;
    mdl_linear.R2adjusted = mdl_linear.fit.Rsquared.Adjusted;
    disp(mdl_linear.R2ordinary)

    % calculate regression line
    % y = intercept + x1*cos(x) + x2 * sin(x)
    x = linspace(-pi, pi, 1000);
    y = mdl_linear.fit.Coefficients.Estimate(1) ...
        + mdl_linear.fit.Coefficients.Estimate(2).* cos(x) ...
        + mdl_linear.fit.Coefficients.Estimate(3).* sin(x);
    mdl_linear.opt = x(y == max(y));
    end

    % Analyse relationship between phase and MEP in moving window approach
    init.T_max = pi;
    init.T_min = -pi;
    current_subject = settings.subjects(idx_participant);
    % generate ground truth for computation grid
    subject_data.phases = phase;
    subject_data.transformed_mep_amplitudes = mep_size';
    simulation_parameters.percentage_data = 0.1;          % Percentage of complete window to be considered
    simulation_parameters.steps =   1000;
    simulation_parameters.wintype = 'winsize';
    simulation_parameters.winsize = (init.T_max - init.T_min)*simulation_parameters.percentage_data;  % window size of simulation
    simulation_parameters.limits =  [init.T_min, init.T_max];
    ground_truth = simulation_ground_truth_slidewindow(simulation_parameters, subject_data, 0, false, false, false);
    ground_truth.opt = ground_truth.X(ground_truth.mean == max(ground_truth.mean));
    clear simulation_parameters subject_data

    % save data in struct
    all_trials.phase        = [all_trials.phase,    phase           ];
    all_trials.MEP_log      = [all_trials.MEP_log,  mep_size        ];
    all_trials.participant  = [all_trials.participant, participant  ];
    all_trials.outlier_all  = [all_trials.outlier_all, participant*0];

    all_results.participant(idx_participant) = settings.subjects(idx_participant);
    all_results.moving_win.opt(idx_participant) = ground_truth.opt(1);
    all_results.c2l_reg.opt(idx_participant) = mdl_linear.opt;
    all_results.c2l_reg.R2ordinary(idx_participant) = mdl_linear.R2ordinary;
    all_results.c2l_reg.p_uncorrected(idx_participant) = mdl_linear.pval;

end % for idx_participant

% Bonferroni correct p-values
all_results.c2l_reg.p_corrected = all_results.c2l_reg.p_uncorrected / length(all_results.participant);

end