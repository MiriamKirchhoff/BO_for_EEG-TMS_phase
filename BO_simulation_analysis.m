%% ANALYSIS SCRIPT FOR BO SIMULATION RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Analysis script for the output files of the Bayesian optimization script
%
% version   1.0, 19.04.2024
% author    Miriam Kirchhoff
% project   C2B
%% Load data

clear
close all
load('\results\simulation_output.mat')


%% Phase location error: based on sliding window and regression

% adapt sliding window result to settings
for i = 1:length(settings.subjects)
    current_subject = settings.subjects(i);
    % generate ground truth for computation grid
    subject_data.phases = all_trials.phase(...
        all_trials.participant == current_subject & ~all_trials.outlier_all);
    subject_data.transformed_mep_amplitudes = all_trials.MEP_log(...
        all_trials.participant == current_subject & ~all_trials.outlier_all)';
    simulation_parameters.percentage_data = 0.1;          % % of complete window to be considered
    simulation_parameters.steps =   n;
    simulation_parameters.wintype = 'winsize';
    simulation_parameters.winsize = (init.T_max - init.T_min)*simulation_parameters.percentage_data;  % window size of simulation
    simulation_parameters.limits =  [init.T_min, init.T_max];
    ground_truth = simulation_ground_truth_slidewindow(simulation_parameters, subject_data, 0, false, false, false);
    [a, b] = max(ground_truth.mean);
    all_results.moving_win.opt(all_results.participant == current_subject) = T(b);
end % for i = 1:length(settings.subjects)
clear a b

% error score: sliding window
analysis_results.phase_error_sliding_window = simulation_results_opt_estimate;
for i = 1:length(settings.subjects)
    % difference to opt
    analysis_results.phase_error_sliding_window(:,:,:,:,i) = ...
        analysis_results.phase_error_sliding_window(:,:,:,:,i) - ...
        all_results.moving_win.opt(all_results.participant == settings.subjects(i));
end

% calculate circular difference
analysis_results.phase_error_sliding_window = ...
    min(abs(analysis_results.phase_error_sliding_window), abs(2*pi - ...
    abs(analysis_results.phase_error_sliding_window)))./pi;

% error score: regression fit
analysis_results.phase_error_regression = simulation_results_opt_estimate;
for i = 1:length(settings.subjects)
    % difference to opt
    analysis_results.phase_error_regression(:,:,:,:,i) = ...
        analysis_results.phase_error_regression(:,:,:,:,i) - ...
        all_results.c2l_reg.opt(all_results.participant == settings.subjects(i));
end % for i = 1:length(settings.subjects)


% calculate circular difference
analysis_results.phase_error_regression = ...
    min(abs(analysis_results.phase_error_regression), abs(2*pi - ...
    analysis_results.phase_error_regression))/pi;


%% MEP size error: Excitability relative difference to optimum

% get X_int from data
simulation_results_opt_estimate_int = ...
    int64((simulation_results_opt_estimate + pi)/(2*pi/199)+1);
analysis_results.excitability_error = NaN(size(simulation_results_opt_estimate_int));
analysis_results.excitability_mean = NaN(size(simulation_results_opt_estimate_int));

% calculate excitability at location
for i = 1:length(settings.subjects)
    current_subject = settings.subjects(i);
    % generate ground truth for computation grid
    subject_data.phases = all_trials.phase(...
        all_trials.participant == current_subject & ~all_trials.outlier_all);
    subject_data.transformed_mep_amplitudes = zscore(all_trials.MEP_log(...
        all_trials.participant == current_subject & ~all_trials.outlier_all))';
    ground_truth = simulation_ground_truth_slidewindow(simulation_parameters, subject_data, 0, false, false, false);
    % calculate the mean excitability at the selected phase
    analysis_results.excitability_mean(:,:,:,:,i) = ...
        ground_truth.mean(simulation_results_opt_estimate_int(:,:,:,:,i));
    % calculate the error score
    % get expected value (= 50% margin)
    ev_gt = mean(ground_truth.mean);
    opt_gt = max(ground_truth.mean);
    % calculate percentage at which the result is dependent on mean being
    % 50% and max being 100%
    percentage_step = (opt_gt - ev_gt)/0.5;
    % (optimum - value) / percentage step is error percentage
    analysis_results.excitability_error(:,:,:,:,i) = ...
        (opt_gt - ground_truth.mean(simulation_results_opt_estimate_int(:,:,:,:,i))) / percentage_step;
end % for i = 1:length(settings.subjects)


%% Calculate convergence

% set how many samples need to not change for convergence criterion
conv_range = 1:20;
criterion = deg2rad(5);

% calculate the earliest convergence
convergence_locations = NaN([length(conv_range) size(simulation_results_opt_estimate, [1 3 4 5])]);

% iterate over the covariance ranges
for idx_range = 1:length(conv_range)
    current_range = conv_range(idx_range);
    % iterate over starting positions
    for start_loc = 1:size(simulation_results_opt_estimate,2)-current_range
        local_change = squeeze(range(simulation_results_opt_estimate(:,start_loc:start_loc+current_range-1,:,:,:), 2));

        new_opt = isnan(squeeze(convergence_locations(idx_range,:,:,:,:))) & local_change < criterion;
        convergence_locations(idx_range,new_opt) = start_loc;

    end % for start_loc = 1:size(simulation_results_opt_estimate,2)

end % for range = 1:length(conv_range)


%% Calculate convergence errors

% init
conv_error.phase = NaN(size(convergence_locations));
conv_error.MEP =  NaN(size(convergence_locations));

% insane for-loop due to easy readability
for i = 1:subject_idx
    disp(i)
    for m = 1:model_idx
        for acq = 1:acquisition_idx
            for reps = 1:length(conv_range)
                for iter = 1:settings.n_repetitions
                    temp_accuracies = analysis_results.phase_error_sliding_window(iter, :, acq, m, i);
                    temp_convergence_locs = convergence_locations(reps, iter, acq, m, i);
                    if ~isnan(temp_convergence_locs)
                        conv_error.phase(reps, iter, acq, m, i) = ...
                            temp_accuracies(temp_convergence_locs);
                    end
                    temp_accuracies = analysis_results.excitability_error(iter, :, acq, m, i);
                    temp_convergence_locs = convergence_locations(reps, iter, acq, m, i);
                    if ~isnan(temp_convergence_locs)
                        conv_error.MEP(reps, iter, acq, m, i) = ...
                            temp_accuracies(temp_convergence_locs);
                    end
                end % for iter = 1:settings.n_repetitions
            end % for reps = 1:length(conv_range)
        end % acq = 1:acquisition_idx
    end % for m = 1:model_idx
end % for i = 1:subject_idx


%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings
clear colors
colors.lightblue    = [45 160 220]/255;
colors.darkgrey     = [1 1 1]*0.4;
colors.lightred     = [220 70 80]/255;
colors.lightgreen   = [220 220 40]/255;
colors.darkred      = [150 30 40]/255;
colors.darkblue     = [40 110 150]/255;
colors.lightgrey    = [1 1 1]*0.5;

colornames = fieldnames(colors);
linestyles = ["--", "-", "-."];
% Define plotting font
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'defaultLineLineWidth',2)

path_current.save = '';
mkdir([path_current.save '\figures\'])

%% Figure 1: Accuracy over trials

% subject_subselection = 'all';
% subject_subselection = {'all', 'median upper', 'median lower'};

accuracy_iteration = figure('units','normalized','outerposition',[0.5 0.5 0.3 .45]);
t = tiledlayout(2, 2, 'TileSpacing','compact');

mdl_names = {'Gaussian Process', 'Bayesian Regression'};
acq_names = {'Adaptive sampling', 'Random sampling'};

% for idx_subject_subselection = 1:numel(subject_subselection)
sub = 1:numel(settings.subjects);
plot_size = [1 1];

iter_end = [100 1000];
ylabels = {'Fast optimization' 'Long-term optimization'};

for idx_iter = 1:numel(iter_end)

    %plot_size = [3-idx_iter 1];

    nexttile(plot_size); grid on; hold on
    acc_current = 1- analysis_results.phase_error_sliding_window(:,:,:,:,sub);

    % calculate parameters
    error_means_SW = squeeze(mean(acc_current, [1 5], 'omitnan'));
    SEM = squeeze(std(acc_current, 0, [1 5], 'omitnan')/sqrt(size(analysis_results.phase_error_sliding_window(:,:,:,:,sub),2)));               % Standard Error
    ts = tinv([0.025  0.975],size(acc_current,2)-1);      % T-Score
    CI_low  = error_means_SW + ts(1)*SEM;
    CI_high = error_means_SW + ts(2)*SEM;                      % two-sided 95% CI
    lgd = [];

    for m = 1:model_idx
        for acq = 1:acquisition_idx
            % plot CI underneath means
            Y = [CI_low(:,acq,m)', flip(CI_high(:,acq,m))'];
            X_area = [1:numel(CI_high(:,acq,m)), flip(1:numel(CI_high(:,acq,m)))];
            fill(X_area, Y, 'm', 'LineStyle', 'none', 'FaceColor', colors.(colornames{acq}), 'FaceAlpha', 0.1)
            lgd = [lgd, ""];
        end % for acq = 1:acquisition_idx
    end % m = 1:model_idx

    for m = 1:model_idx
        for acq = 1:acquisition_idx
            % plot mean
            plot(error_means_SW(:,acq,m), 'Color', colors.(colornames{acq}),...
                'linewidth', 2, 'linestyle', linestyles(m))
            lgd = [lgd, [settings.type_acquisition{acq} ', ' settings.type_model{m}]];
        end % for acq = 1:acquisition_idx
    end % m = 1:model_idx

    xlim([1, iter_end(idx_iter)])%settings.n_iterations+1])
    if idx_iter == 1; title('Phase location accuracy ', 'FontName', 'Times New Roman'); end
    ylabel(ylabels{idx_iter})

    nexttile(plot_size); grid on; hold on

    acc_current = 1 - analysis_results.excitability_error(:,:,:,:,sub);

    % calculate parameters
    error_means_SW = squeeze(mean(acc_current, [1 5], 'omitnan'));
    SEM = squeeze(std(acc_current, 0, [1 5], 'omitnan')/sqrt(size(analysis_results.excitability_error(:,:,:,:,sub),2)));               % Standard Error
    ts = tinv([0.025  0.975],size(acc_current,2)-1);      % T-Score
    CI_low  = error_means_SW + ts(1)*SEM;
    CI_high = error_means_SW + ts(2)*SEM;                      % two-sided 95% CI
    lgd = [];

    for m = 1:model_idx
        for acq = 1:acquisition_idx
            % plot CI underneath means
            Y = [CI_low(:,acq,m)', flip(CI_high(:,acq,m))'];
            X_area = [1:numel(CI_high(:,acq,m)), flip(1:numel(CI_high(:,acq,m)))];
            fill(X_area, Y, 'm', 'LineStyle', 'none', 'FaceColor', colors.(colornames{acq}), 'FaceAlpha', 0.1)
            lgd = [lgd, ""];
        end % for acq = 1:acquisition_idx
    end % m = 1:model_idx

    for m = flip(1:model_idx) % flip for legend order
        plot(NaN,NaN,'w');
        lgd = [lgd, ['\textbf{' mdl_names{m} '}']];
        for acq = 1:acquisition_idx
            % plot mean
            plot(error_means_SW(:,acq,m), 'Color', colors.(colornames{acq}),...
                'linewidth', 2, 'linestyle', linestyles(m))
            lgd = [lgd, [acq_names{acq}]];
        end % for acq = 1:acquisition_idx
    end % m = 1:model_idx

    xlim([1, iter_end(idx_iter)])%settings.n_iterations+1])
    if idx_iter == 1; title('MEP size accuracy ', 'FontName', 'Times New Roman'); end

end

xlabel(t, 'Iteration', 'FontName', 'Times New Roman')
ylabel(t, 'Accuracy', 'FontName', 'Times New Roman')

title(t, 'Progression of accuracy over iterations ', 'FontName', 'Times New Roman')
l = legend(lgd, 'Location', 'southoutside', 'NumColumns', 2, 'Orientation','Vertical', 'interpreter','latex');
l.Layout.Tile = 'north';
%fontsize(gcf,scale=1.4)
saveas(accuracy_iteration,[path_current.save '\figures\3iteration.jpg'])


%% Display accuracy at tipping points

% define iterations to be calculated
iter = [1 101 1001];

fprintf('\n \n Phase accuracy \n \n')
% phase error
for idx_iter = 1:numel(iter)
    if length(iter) <  iter(idx_iter)
        disp(['calculating iteration accuracy not possible.' ...
            'simulation contains less iterations than criterion'])
        continue
    end
    for idx_mdl = [2 1]
        for idx_acq = 1:acquisition_idx
            acc_current = 1 - analysis_results.phase_error_sliding_window(:,iter(idx_iter),idx_acq,idx_mdl,sub);

            error_means_SW = squeeze(mean(acc_current, [1 5], 'omitnan'));
            SEM = squeeze(std(acc_current, 0, [1 5], 'omitnan')/sqrt(size(analysis_results.phase_error_sliding_window(:,:,:,:,sub),2)));               % Standard Error
            ts = tinv([0.025  0.975],length(acc_current)-1);      % T-Score
            CI_low  = error_means_SW + ts(1)*SEM;
            CI_high = error_means_SW + ts(2)*SEM;                      % two-sided 95% CI

            fprintf('iter = %d, %s, %s \n', iter(idx_iter), acq_names{idx_acq}, mdl_names{idx_mdl})
            fprintf('%.1f  [%.1f %.1f] \n', error_means_SW*100, CI_low*100, CI_high*100)

        end % for idx_acq = 1:acquisition_idx
    end % for idx_mdl = 1:model_idx
end % for idx_iter = 1:numel(iter)

fprintf('\n \n MEP accuracy \n \n')

% phase error
for idx_iter = 1:numel(iter)
    if length(iter) <  iter(idx_iter)
        disp(['calculating iteration accuracy not possible.' ...
            'simulation contains less iterations than criterion'])
        continue
    end
    for idx_mdl = [2 1]
        for idx_acq = 1:acquisition_idx
            acc_current = 1 - analysis_results.excitability_error(:,iter(idx_iter),idx_acq,idx_mdl,sub);

            error_means_SW = squeeze(mean(acc_current, [1 5], 'omitnan'));
            SEM = squeeze(std(acc_current, 0, [1 5], 'omitnan')/sqrt(size(analysis_results.phase_error_sliding_window(:,:,:,:,sub),2)));               % Standard Error
            ts = tinv([0.025  0.975],length(acc_current)-1);      % T-Score
            CI_low  = error_means_SW + ts(1)*SEM;
            CI_high = error_means_SW + ts(2)*SEM;                      % two-sided 95% CI

            fprintf('iter = %d, %s, %s \n', iter(idx_iter), acq_names{idx_acq}, mdl_names{idx_mdl})
            fprintf('%.1f  [%.1f %.1f] \n', error_means_SW*100, CI_low*100, CI_high*100)

        end % for idx_acq = 1:acquisition_idx
    end % for idx_mdl = 1:model_idx
end % for idx_iter = 1:numel(iter)



%% Figure 2: Convergence criterion selection

sub = 1:numel(settings.subjects);
plot_size = [1 1];

conv_plot = figure('units','normalized','outerposition',[0.5 0.5 0.3 .4]);
t = tiledlayout('flow', 'TileSpacing', 'compact');
title(t, 'Convergence analysis', 'FontName', 'Times New Roman')

error_type = fieldnames(conv_error);
error_name = {'Phase location', 'MEP size'};

% Subplot: Error vs nr consecutive trials
for idx_error_type = 1:numel(error_type)
    nexttile; grid on; hold on
    title([error_name{idx_error_type} ' accuracy'])
    for m = [2 1]
        for acq = 1:acquisition_idx
            plot(NaN,NaN,'w')
            acc_current = 1 - conv_error.(error_type{idx_error_type})(:,:,acq, m, sub);
            temp_mean = mean(acc_current, [2 5], 'omitnan');
            % calculate CI
            N = size(acc_current,2)*size(conv_error.(error_type{idx_error_type})(:,:,acq, m, sub),5);
            SEM = squeeze(std(acc_current, 0, [2 5], "omitnan")/sqrt(N));
            ts = tinv([0.025  0.975], N-1);      % T-Score
            CI_low  = temp_mean + ts(1)*SEM;
            CI_high = temp_mean + ts(2)*SEM;

            % plot 95% CI
            Y = [CI_low', flip(CI_high)'];
            X_area = [1:numel(CI_high), flip(1:numel(CI_high))];
            fill(X_area, Y, 'm', 'LineStyle', 'none', 'FaceColor', colors.(colornames{acq}), 'FaceAlpha', 0.2)
        end
        for acq = 1:acquisition_idx
            acc_current = 1 - conv_error.(error_type{idx_error_type})(:,:,acq, m, sub);
            temp_mean = mean(acc_current, [2 5], 'omitnan');
            % plot mean
            plot(temp_mean, 'Color', colors.(colornames{acq}), 'linewidth', 2, ...
                'linestyle', linestyles(m))
        end
    end % m = 1:model_idx


    if idx_error_type == 1; ylabel('Accuracy'); end
    % ylim([0.6 0.85])
    xlim([1 20])
end

% Subplot: Convergence iteration vs criterion
nexttile([1 2]); grid on; hold on
title('Iteration of convergence')
hold on
grid on
for m = [2 1]
    plot(NaN,NaN,'w')
    for acq = 1:acquisition_idx

        temp_mean = mean(convergence_locations(:,:,acq, m, sub), [2 5], 'omitnan');
        % calculate CI
        N = size(convergence_locations(:,:,acq, m, sub),2)*size(convergence_locations(:,:,acq, m, sub),5);
        SEM = squeeze(std(convergence_locations(:,:,acq, m, sub), 0, [2 5], "omitnan")/sqrt(N));
        ts = tinv([0.025  0.975], N-1);      % T-Score
        CI_low  = temp_mean + ts(1)*SEM;
        CI_high = temp_mean + ts(2)*SEM;

        % plot 95% CI
        Y = [CI_low', flip(CI_high)'];
        X_area = [1:numel(CI_high), flip(1:numel(CI_high))];
        fill(X_area, Y, 'm', 'LineStyle', 'none', 'FaceColor', colors.(colornames{acq}), 'FaceAlpha', 0.2)
    end
    for acq = 1:acquisition_idx
        temp_mean = mean(convergence_locations(:,:,acq, m, sub), [2 5], 'omitnan');
        % plot mean
        plot(temp_mean, 'Color', colors.(colornames{acq}), 'linewidth', 2, ...
            'linestyle', linestyles(m))
    end
end % m = 1:model_idx

ylabel('Iteration')
xlim([min(conv_range)-0.5, max(conv_range)+0.5])
xlim([1 20])

xlabel(t, 'Convergence criterion', 'FontName', 'Times New Roman')
lgd = ["\textbf{Bayesian Regression}" "" "" "Adaptive sampling" "Random sampling" "\textbf{Gaussian Process}" "" "" "Adaptive sampling" "Random sampling"];
l = legend(lgd, 'Location', 'southoutside', 'NumColumns', 2, 'Orientation','Vertical', 'interpreter','latex');
l.Layout.Tile = 'north';
