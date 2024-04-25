function [ground_truth, shuffle] = simulation_ground_truth_slidewindow(simulation_parameters, subject_data, shuffle_repetitions, robust_estimate, correct_CI, plot_figs)
% function [ground_truth, shuffle] = simulation_ground_truth_slidewindow(simulation_parameters,subject_data, shuffle_repetitions)
%
% Calculates for entire array with simulation_parameter.steps steps the ground truth in terms
% of mean, standard deviation, and confidence interval. Plots the data,
% mean, std.
%
% INPUTS
% simulation_parameters: Struct with fields:
    % steps:    Number of samples being evaluated, integer
    % wintype:  define whether the window is limited by number of samples
    %           ('sample_absolute'), percentage of total samples in decimals
    %           ('sample_relative'), or total windowsize ('winsize'), string
    % winsize:  define window size based on windowtype
    % limits:   define min and max period of the data as vector [min max]
% subject_data: struct with fields:
    % phases:   phase raw data, vector 1xn
    % transformed_mep_amplitudes: mep raw data, vector 1xn
% shuffle_repetitions: number of repetitions for generation of
%               surrogate data. No creation of surrogate data if
%               shuffle_repetitions = 0. Integer.
% robust_estimate: If true, replace mean and std by median and mad,
%           default is false, bool.
% correct_CI: If no entry, this is true. controls whether or not
%               to apply multiple testing correction to the CI. Bool.
%
% OUTPUTS
% ground_truth: Struct with fields:
    % X:        x-values at which the data was evaluated using the
    %           sliding window,         vector 1xlength(X)
    % mean:     mean of the data at X,  vector 1xlength(X)
    % sd:       sd of the data at X,    vector 1xlength(X)
    % ci:       ci of the data at X,    vector 1xlength(X)
    % shuffle:      Struct with fields:
    % mean:     mean of the shuffeled data, double
    % ci:       ci of the shuffeled data, double
%
% version   1.0, 21.02.2023
% author    Miriam Kirchhoff
% project   C2B

%%
if nargin < 4
    robust_estimate = false;
end

if nargin < 5
    correct_CI = true;
end

if nargin < 6
    plot_figs = true;
end

simulation_parameters.dataset = [subject_data.phases',subject_data.transformed_mep_amplitudes];


%% Calculate sliding window over all ROI

% determine all locations to evaluate
x_values = linspace(min(simulation_parameters.limits),max(simulation_parameters.limits),simulation_parameters.steps);

% initialize
mean_simulated_mep  = nan(1,simulation_parameters.steps);
std_simulated_mep   = nan(1,simulation_parameters.steps);
ci_simulated_mep    = nan(2,simulation_parameters.steps);

% iterate over locations to evaluate
for i = 1:simulation_parameters.steps
    [~, mean_simulated_mep(i), std_simulated_mep(i), ci_simulated_mep(:,i)] = ...
        simulate_data_slidewindow(x_values(i), simulation_parameters.dataset, simulation_parameters.wintype, simulation_parameters.winsize, simulation_parameters.limits, robust_estimate);
end % i = 1:simulation_parameters.steps

% save output variables
ground_truth.X = x_values;
ground_truth.mean = mean_simulated_mep;
ground_truth.sd = std_simulated_mep;
ground_truth.ci = ci_simulated_mep;


%% Create surrogate data

% disp('start generation of surrogate data')

% initialize
mean_shuffle_mep = nan(simulation_parameters.steps, shuffle_repetitions);
ci_shuffle_mep_lower = nan(simulation_parameters.steps, shuffle_repetitions);
ci_shuffle_mep_higher = nan(simulation_parameters.steps, shuffle_repetitions);

for shuffle_iteration = 1:shuffle_repetitions   % loop for number of shuffles required

    % shuffle data
    dataset_shuffle = [subject_data.phases',subject_data.transformed_mep_amplitudes(randperm(length(subject_data.transformed_mep_amplitudes)))];

    for i = 1:simulation_parameters.steps % loop over locations of interest
        if correct_CI
            [~, mean_shuffle_mep(i,shuffle_iteration), ~, ci_shuffle_mep] = simulate_data_slidewindow(x_values(i), dataset_shuffle, simulation_parameters.wintype, simulation_parameters.winsize, simulation_parameters.limits, robust_estimate, simulation_parameters.steps);
        else
            [~, mean_shuffle_mep(i,shuffle_iteration), ~, ci_shuffle_mep] = simulate_data_slidewindow(x_values(i), dataset_shuffle, simulation_parameters.wintype, simulation_parameters.winsize, simulation_parameters.limits, robust_estimate);
        end
        ci_shuffle_mep_lower(i,shuffle_iteration) = min(ci_shuffle_mep);
        ci_shuffle_mep_higher(i,shuffle_iteration) = max(ci_shuffle_mep);
    end % i = 1:simulation_parameters.steps
end % for shuffle_iteration = 1:shuffle_repetitions

% Calculate mean data across shuffle iterations and locations
ci_shuffle_mep_higher = mean(ci_shuffle_mep_higher,'all');
ci_shuffle_mep_lower = mean(ci_shuffle_mep_lower,'all');

% save output variables for shuffle
shuffle.ci = [ci_shuffle_mep_higher, ci_shuffle_mep_lower];
shuffle.mean = mean(mean_shuffle_mep,'all');    % should be the overall mean of the data

    	
if plot_figs
    %% Plot results

    % t = tiledlayout(1,1);
    % nexttile
    grid on

    % plot all datapoints
    scatter(subject_data.phases, subject_data.transformed_mep_amplitudes, 'k', 'filled','o','MarkerEdgeColor','none','MarkerFaceAlpha',0.1)
    hold on

    % plot mean line
    plot(x_values,mean_simulated_mep,'k','LineWidth',2)

    % plot standard deviation
    X_area = [x_values, flip(x_values)];
    Y = [mean_simulated_mep+std_simulated_mep, flip(mean_simulated_mep-std_simulated_mep)];
    fill(X_area, Y, 'm', 'LineStyle', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2)

    % plot oscillation w.r.t. phase
    plot(x_values, mean(std_simulated_mep)*cos(x_values) + min(subject_data.transformed_mep_amplitudes)+1, ':k', LineWidth=2)
    max_amp = max(subject_data.transformed_mep_amplitudes);
    max_mep = min(x_values(mean_simulated_mep == max(mean_simulated_mep))); % use min in case max is at +/- pi

    % plot windowsize
    line([max_mep-0.5*simulation_parameters.winsize, max_mep-0.5*simulation_parameters.winsize, max_mep+0.5*simulation_parameters.winsize, max_mep+0.5*simulation_parameters.winsize,], ...
        [max_amp, max_amp+0.1*mean(std_simulated_mep), max_amp+0.1*mean(std_simulated_mep), max_amp],...
        'Color', [0.5 0.5 0.5], 'LineWidth', 2)

    % plot max mep location
    xline(max_mep,'-',{'Maximum MEP   '})

    xlabel('phase [rad]')
    ylabel('MEP amplitude (z-score)')
    xlim(simulation_parameters.limits)

    % plot surrogate data analytical
    Y = [ones(size(x_values)).*ci_shuffle_mep_higher, flip(ones(size(x_values)).*ci_shuffle_mep_lower)];
    fill(X_area, Y, 'r', 'LineStyle', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.2)

    alpha = 0.05;
    if correct_CI
        alpha = alpha/simulation_parameters.steps;
    end

end % if plot

end % end of function