%% Simulation parameter extraction for preprocessed data
%
% version   1.0, 19.04.2024
% author    Miriam Kirchhoff
% project   C2B

clear
close all
clc


%% Settings

subjects = 1:38;

path_current.load_EEG = '\data_preprocessed';
path_current.load_EMG = '\data_epoched';
path_current.save = '\simulation_data';
mkdir(path_current.save)
formatSpec = '%03.0f';

% settings for montage
hjorth.channel = 'C3';

% setttings for EMG montage
merge_method = 'pca';


%% Initialize
all_results = [];
all_trials = [];
initial_subject = true;
plotting = true;


%% Iterate over subjects

for idx_subject = subjects  % iterate over subjects

    clearvars -except idx_subject phast all_results all_trials ...
        plotting initial_subject path_current formatSpec hjorth merge_method


    %% Load preprocessed EEG data

    % start timer
    tic

    % start notification
    fprintf('\nStarting participant %03.0f. \n', idx_subject)
    disp('Start loading...')

    % load EEG data
    load([path_current.load_EEG '\' num2str(idx_subject,formatSpec) '\EEG_lap_P' num2str(idx_subject,formatSpec)]);
    
    % timing notification
    fprintf('loading completed. Expired time: %.0f seconds \n', toc)


    %% Calculate bandpower

    % Can specify band here if of special interest (e.g. mu)
    parameters.bandpower = bandpower(EEG_lap.data);


    %% Calculate phase

    % Phastimate settings
    phast.edge = 65;
    phast.ord = 30;
    phast.hilbertwindow = 128;
    phast.offset_correction = 4;

    parameters.phase = phastimate_adapted(EEG_lap.data, phast.edge, ...
        phast.ord, phast.hilbertwindow, phast.offset_correction);


    %% Calculate MEP

    % settings: time windows for preinnervation and MEP
    mep_amp.prestim_time = [-100 -5];
    mep_amp.limits = [20 40];

    % load EMG data
    load([path_current.load_EMG '\' num2str(idx_subject,formatSpec) '\EMG_P' num2str(idx_subject,formatSpec)]);
    EMG = EMG_epoched;
    clear EMG_epoched

    % baseline removal for plotting
    EMG = mk_rmbase(EMG, mep_amp.limits);

    % remove EEG outlier trials
    EMG = pop_select(EMG, 'trial', find(~rejection.threshold));

    % detrending
    EMG_detrended = tesa_detrend(EMG, 'linear', mep_amp.limits);

    % merge data of two channels
    switch merge_method % how to merge EMG channels
        case 'max'      % merge based on higher value
            emg_merged.raw = squeeze(max(EMG.data));
            emg_merged.detrended = squeeze(max(EMG_detrended.data));
        case 'pca'      % merge based on Zrenner (2023) PCA method
            % calculate coefficients of pca
            coeff = pca(reshape(EMG.data, size(EMG.data,1), [])');
            % map onto first coefficient
            temp = reshape(EMG.data, size(EMG.data,1), [])'*coeff(:,1);
            % reshape back into trials
            emg_merged.raw = reshape(temp, size(EMG.data,2), [])';

            % calculate coefficients of pca
            coeff = pca(reshape(EMG_detrended.data, size(EMG_detrended.data,1), [])');
            % map onto first coefficient
            temp = reshape(EMG_detrended.data, size(EMG_detrended.data,1), [])'*coeff(:,1);
            % reshape back into trials
            emg_merged.detrended = reshape(temp, size(EMG_detrended.data,2), [])';
    end

    mep_data = emg_merged.raw(:, EMG.times >= mep_amp.limits(1) & ...
        EMG.times <= mep_amp.limits(2), :);
    pi_data = emg_merged.detrended(:, EMG.times >= mep_amp.prestim_time(1) & ...
        EMG.times <= mep_amp.prestim_time(2),:);


    %% Calculate MEPs: find all peaks in data and get range of that

    prominence = 1;

    for i = 1:size(mep_data,1)
        % find maximum peaks
        [pks_max, temp_loc_max] = findpeaks(mep_data (i,:), 'MinPeakProminence',prominence);
        % find minimum peaks
        [pks_min, temp_loc_min]  = findpeaks(-mep_data (i,:), 'MinPeakProminence',prominence);

        % if at least one peak was not found, enter NaN
        if isempty(pks_min) || isempty(pks_max) % if one peak was not found
            parameters.mep_raw(i) = NaN;
            loc_max(i) = NaN;
            loc_min(i) = NaN;
        else
            % get maximum range and location of the highest peaks
            loc_max(i) = temp_loc_max(find(pks_max == max(pks_max), 1, 'first'));
            loc_min(i) = temp_loc_min(find(pks_min == max(pks_min),1,'first'));
            parameters.mep_raw(i) = max(pks_max) - min(-pks_min);
        end
    end

    % log scale
    parameters.mep_log = log(parameters.mep_raw)';

    % z-score data
    parameters.mep_log = parameters.mep_log - mean(parameters.mep_log, 'omitnan');
    parameters.mep_log = parameters.mep_log / std(parameters.mep_log, 'omitnan');

    parameters.preinnervation = squeeze(range(pi_data,2))';


    %% Outlier rejection EMG

    % Preinnervation cutoff
    parameters.pi_cutoff = 50;
    rejection.preinnervation = (parameters.preinnervation > parameters.pi_cutoff);
    rejection.no_peak = isnan(parameters.mep_raw) | parameters.mep_raw < parameters.pi_cutoff;
    rejection.all_mep_criteria = rejection.preinnervation | rejection.no_peak;


    %% Investigate local trends

    slidingmedian.size = 100;       % number of samples to be taken into account

    % sliding window median
    for i = 1:length(parameters.mep_raw)
        min_loc = max(1, i-round(0.5*slidingmedian.size)-1);
        max_loc = min(length(parameters.mep_raw), i+round(0.5*slidingmedian.size));
        temp = min_loc:max_loc;
        temp(rejection.all_mep_criteria(temp)) = [];
        slidingmedian.median(i) = median(parameters.mep_raw(temp), "omitnan");
        slidingmedian.mad(i) = mad(parameters.mep_raw(temp));
        if slidingmedian.mad(i) == 0 % for snippets with only outliers
            slidingmedian.mad(i) = slidingmedian.mad(i-1);
        end
    end

    parameters.mep_slidemed = parameters.mep_raw-slidingmedian.median;
    % parameters.mep_slidemed = parameters.mep_slidemed./slidingmedian.mad;

    % t-score data (mean 50, sd 10) so there are likely no negative numbers
    % for log transform
    parameters.mep_slidemed = (parameters.mep_slidemed / std(parameters.mep_slidemed, 'omitnan'));
    parameters.mep_slidemed = parameters.mep_slidemed - min(parameters.mep_slidemed) + .5;

    % log transform data
    parameters.mep_slidemed = log(parameters.mep_slidemed);

    % z-score results
    parameters.mep_slidemed = (parameters.mep_slidemed / std(parameters.mep_slidemed, 'omitnan')) - mean(parameters.mep_slidemed, 'omitnan');


    %% Analysis

    [rho, pval] = circ_corrcl(parameters.phase(~rejection.all_mep_criteria), parameters.mep_slidemed(~rejection.all_mep_criteria));

    mdl_linear.estphase = parameters.phase(~rejection.all_mep_criteria);    % enter phase values
    mdl_linear.excitability = zscore(parameters.mep_slidemed(~rejection.all_mep_criteria));     % enter MEPs
    mdl_linear.fit = fitlm([cos(mdl_linear.estphase)' sin(mdl_linear.estphase)'], mdl_linear.excitability);

    mdl_linear.pval = mdl_linear.fit.coefTest;
    mdl_linear.R2ordinary = mdl_linear.fit.Rsquared.Ordinary;
    mdl_linear.R2adjusted = mdl_linear.fit.Rsquared.Adjusted;

    % calculate regression line
    % y = intercept + x1*cos(x) + x2 * sin(x)
    x = linspace(-pi, pi, 1000);
    y = mdl_linear.fit.Coefficients.Estimate(1) ...
        + mdl_linear.fit.Coefficients.Estimate(2).* cos(x) ...
        + mdl_linear.fit.Coefficients.Estimate(3).* sin(x);
    mdl_linear.opt = x(y == max(y));

    mg_step = 0.05;                 % measurement grid stepsize
    min_phase = -pi;                % min possible phase
    max_phase = pi;                 % max possible phase
    subject_data.transformed_mep_amplitudes = zscore(parameters.mep_slidemed(~rejection.all_mep_criteria))';
    subject_data.phases = parameters.phase(~rejection.all_mep_criteria);
    simulation_parameters.percentage_data = 0.2;
    simulation_parameters.steps =   100;
    simulation_parameters.dataset = [subject_data.phases',subject_data.transformed_mep_amplitudes];
    simulation_parameters.wintype = 'winsize';
    simulation_parameters.winsize = (max_phase-min_phase)*simulation_parameters.percentage_data;  % window size of simulation
    simulation_parameters.limits =  [min_phase, max_phase];
    shuffle_repetitions =           5000;
    cluster_based_perm_sig = cluster_perm_test_slidingwindow(simulation_parameters,subject_data, shuffle_repetitions, 0.05/32);


    %% Plotting EMG information

    if plotting

        trials_plot = figure('units','normalized','outerposition',[0.5 0 0.5 1]);

        % plot rejection histogram
        nexttile; hold on; grid on
        h = histogram('binedges', [0:3], 'bincounts', [sum(rejection.threshold) sum(rejection.preinnervation) sum(rejection.no_peak)], 'FaceColor','k');
        % write numbers on top
        y = h.BinCounts;
        x = 0.5:1:2.5;
        text(x, y+50, string(y))
        xticks(x)
        xticklabels({'EEG threshold','preinnervation','no peak'})
        ylim([0 1200])
        ylabel('number of trials rejected')
        title(['Rejection across subject ' num2str(idx_subject)])

        % plot bandpower
        nexttile
        x = 1:length(parameters.bandpower);
        semilogy(x(~rejection.all_mep_criteria), parameters.bandpower(~rejection.all_mep_criteria), 'k.')
        hold on; grid on
        semilogy(x(rejection.preinnervation), parameters.bandpower(rejection.preinnervation), 'c.')
        semilogy(x(rejection.no_peak), parameters.bandpower(rejection.no_peak), 'm.')
        yline(median(parameters.bandpower), 'b', 'median bandpower', 'LineWidth',2)
        ylabel('bandpower')
        xlabel('trial number')
        title(['Distribution of bandpower across subject ' num2str(idx_subject)])

        % plot meps
        nexttile; hold on; grid on
        x = 1:length(parameters.mep_slidemed);
        plot(x(~rejection.all_mep_criteria), parameters.mep_slidemed(~rejection.all_mep_criteria), 'k.', 'linewidth', 1)
        plot(x(rejection.preinnervation), parameters.mep_slidemed(rejection.preinnervation), 'c.', 'linewidth', 1)
        plot(x(rejection.no_peak), parameters.mep_slidemed(rejection.no_peak), 'm.', 'linewidth', 1)
        %yline(parameters.pi_cutoff, 'm', 'outlier threshold', 'LineWidth',2);
        ylabel('MEP (z-score)')
        xlabel('trial number')
        title(['Distribution of MEP across subject ' num2str(idx_subject)])

        % plot preinnervation
        nexttile
        x = 1:length(parameters.preinnervation);
        semilogy(x(~rejection.all_mep_criteria), parameters.preinnervation(~rejection.all_mep_criteria), 'k.')
        hold on
        semilogy(x(rejection.preinnervation), parameters.preinnervation(rejection.preinnervation), 'c.')
        semilogy(x(rejection.no_peak), parameters.preinnervation(rejection.no_peak), 'm.')
        grid on
        yline(parameters.pi_cutoff, 'c', 'outlier threshold', 'LineWidth',2);
        ylabel('preinnervation [µV]')
        xlabel('trial number')
        title(['Distribution of preinnervation across subject ' num2str(idx_subject)])

        % plot EMG response in peak-to-peak window
        nexttile([2 1]); hold on; grid on;
        x = EMG.times(EMG.times >= mep_amp.limits(1) & ...
            EMG.times <= mep_amp.limits(2));
        y = mep_data';
        plot(x, y(:, ~rejection.all_mep_criteria), 'Color', [0 0 0 0.2])
        try
            plot(x, y(:, rejection.preinnervation), 'Color', [0 1 1 0.2])
        end
        try %#ok<TRYNC>
            plot(x, y(:, rejection.no_peak), 'Color', [1 0 1 0.2])
        end
        ylabel('EMG response [µV]')
        xlabel('t w.r.t. TMS pulse [ms]')
        title(['EMG response for subject ' num2str(idx_subject)])

        % show peak locations
        temp1 = sub2ind(size(y),loc_max(~rejection.all_mep_criteria), find(~rejection.all_mep_criteria));
        scatter(x(loc_max(~rejection.all_mep_criteria)), y(temp1), '^k', 'filled', 'SizeData', 15)
        temp2 = sub2ind(size(y),loc_min(~rejection.all_mep_criteria), find(~rejection.all_mep_criteria));
        scatter(x(loc_min(~rejection.all_mep_criteria)), y(temp2), 'vk', 'filled', 'SizeData', 15)
        ylim([min(y(temp2))+min(y(temp2))/30 max(y(temp1))+max(y(temp1))/30])

        % plot simulation data
        nexttile([2 1]); hold on; grid on;
        mg_step = 0.05;                 % measurement grid stepsize
        min_phase = -pi;                % min possible phase
        max_phase = pi;                 % max possible phase
        subject_data.transformed_mep_amplitudes = zscore(parameters.mep_slidemed(~rejection.all_mep_criteria))';
        subject_data.phases = parameters.phase(~rejection.all_mep_criteria);
        simulation_parameters.percentage_data = 0.2;
        simulation_parameters.steps =   100;
        simulation_parameters.dataset = [subject_data.phases',subject_data.transformed_mep_amplitudes];
        simulation_parameters.wintype = 'winsize';
        simulation_parameters.winsize = (max_phase-min_phase)*simulation_parameters.percentage_data;  % window size of simulation
        simulation_parameters.limits =  [min_phase, max_phase];
        shuffle_repetitions =           5;
        [ground_truth, shuffle] = simulation_ground_truth_slidewindow(simulation_parameters, subject_data, shuffle_repetitions, 0, false);

        title(['MEP (no outliers) per phase for subject ' num2str(idx_subject)])
        subtitle(['Rho = ' num2str(rho) ', p = ' num2str(pval)])
        xlim([-pi pi])

        % plot regression line
        % y = intercept + x1*cos(x) + x2 * sin(x)
        x = linspace(-pi, pi, 1000);
        y = mdl_linear.fit.Coefficients.Estimate(1) ...
            + mdl_linear.fit.Coefficients.Estimate(2).* cos(x) ...
            + mdl_linear.fit.Coefficients.Estimate(3).* sin(x);

        plot(x,y, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2)

        legend(["" "sliding window fit +-SD" "" "corresponding oscillation"...
            "" "" "95% CI" "regression fit"])


        %% Save Plotting

        saveas(trials_plot,[path_current.save '\figures\singlesubject_' num2str(idx_subject,formatSpec) '.jpg'])
        close all

    else % if not plotting
        mg_step = 0.05;                 % measurement grid stepsize
        min_phase = -pi;                % min possible phase
        max_phase = pi;                 % max possible phase
        subject_data.transformed_mep_amplitudes = zscore(parameters.mep_slidemed(~rejection.all_mep_criteria))';
        subject_data.phases = parameters.phase(~rejection.all_mep_criteria);
        simulation_parameters.percentage_data = 0.2;
        simulation_parameters.steps =   100;
        simulation_parameters.dataset = [subject_data.phases',subject_data.transformed_mep_amplitudes];
        simulation_parameters.wintype = 'winsize';
        simulation_parameters.winsize = (max_phase-min_phase)*simulation_parameters.percentage_data;  % window size of simulation
        simulation_parameters.limits =  [min_phase, max_phase];
        shuffle_repetitions =           5;
        [ground_truth, shuffle] = simulation_ground_truth_slidewindow(simulation_parameters, subject_data, shuffle_repetitions, 0, false, false);
    end % if plotting


    %% Merge variables

    if initial_subject

        % initialize structure

        % trial-wise data
        all_trials.participant      = ones(size(parameters.bandpower))*idx_subject;
        all_trials.trial            = 1:length(parameters.bandpower);
        all_trials.phase            = parameters.phase;
        all_trials.bandpower        = parameters.bandpower;
        all_trials.MEP_raw          = parameters.mep_raw;
        all_trials.MEP_log          = parameters.mep_log';
        all_trials.preinnervation   = parameters.preinnervation;
        all_trials.outlier_preinnervation = parameters.preinnervation > parameters.pi_cutoff;
        all_trials.outlier_nopeak   = parameters.mep_raw < parameters.pi_cutoff | isnan(parameters.mep_raw);
        all_trials.outlier_all      = all_trials.outlier_preinnervation | all_trials.outlier_nopeak;

        % subject-wise data
        all_results.participant         = idx_subject;
        all_results.trials              = length(parameters.bandpower);
        all_results.normality_violated  = mdl_linear.normality_violated;

        % correlation
        all_results.c2lcorr.p_uncorrected = pval;
        all_results.c2lcorr.p_corrected = NaN;
        all_results.c2lcorr.rho         = rho;

        % rejection
        all_results.rejection.total     = sum(all_trials.outlier_all(all_trials.participant == idx_subject)) + sum(rejection.threshold);
        all_results.rejection.EEG_threshold = sum(rejection.threshold);
        all_results.rejection.preinnervation = sum(rejection.preinnervation);
        all_results.rejection.no_peak   = sum(rejection.no_peak);

        all_results.moving_win.opt      = ground_truth.X(find(ground_truth.mean == max(ground_truth.mean), 1, 'first'));

        all_results.c2l_reg.p_uncorrected = mdl_linear.pval;
        all_results.c2l_reg.p_corrected = NaN;
        all_results.c2l_reg.R2ordinary  = mdl_linear.R2ordinary;
        all_results.c2l_reg.R2adjusted  = mdl_linear.R2adjusted;
        all_results.c2l_reg.intercept   = mdl_linear.fit.Coefficients.Estimate(1);
        all_results.c2l_reg.x1          = mdl_linear.fit.Coefficients.Estimate(2);
        all_results.c2l_reg.x2          = mdl_linear.fit.Coefficients.Estimate(3);
        all_results.c2l_reg.opt         = mdl_linear.opt;

        all_results.cluster_based_perm_sig = cluster_based_perm_sig;

        initial_subject = false;

    else % add to existing data

        all_trials.participant          = [all_trials.participant, ones(size(parameters.bandpower))*idx_subject];
        all_trials.trial                = [all_trials.trial, 1:length(parameters.bandpower)];
        all_trials.phase                = [all_trials.phase, parameters.phase];
        all_trials.bandpower            = [all_trials.bandpower, parameters.bandpower];
        all_trials.MEP_raw              = [all_trials.MEP_raw, parameters.mep_raw];
        all_trials.MEP_log              = [all_trials.MEP_log, parameters.mep_log'];
        all_trials.preinnervation       = [all_trials.preinnervation, parameters.preinnervation];
        all_trials.outlier_preinnervation = [all_trials.outlier_preinnervation, parameters.preinnervation > parameters.pi_cutoff];
        all_trials.outlier_nopeak       = [all_trials.outlier_nopeak, parameters.mep_raw < parameters.pi_cutoff | isnan(parameters.mep_raw)];
        all_trials.outlier_all          = all_trials.outlier_preinnervation | all_trials.outlier_nopeak;

        all_results.participant         = [all_results.participant, idx_subject];
        all_results.trials              = [all_results.trials, length(parameters.bandpower)];

        % correlation
        all_results.c2lcorr.p_uncorrected = [all_results.c2lcorr.p_uncorrected, pval];
        all_results.c2lcorr.p_corrected = [all_results.c2lcorr.p_corrected, NaN];
        all_results.c2lcorr.rho         = [all_results.c2lcorr.rho, rho];

        % rejection
        all_results.rejection.total     = [all_results.rejection.total, sum(all_trials.outlier_all(all_trials.participant == idx_subject)) + sum(rejection.threshold)];
        all_results.rejection.EEG_threshold = [all_results.rejection.EEG_threshold, sum(rejection.threshold)];
        all_results.rejection.preinnervation = [all_results.rejection.preinnervation, sum(rejection.preinnervation)];
        all_results.rejection.no_peak   = [all_results.rejection.no_peak, sum(rejection.no_peak)];

        all_results.moving_win.opt      = [all_results.moving_win.opt, ground_truth.X(find(ground_truth.mean == max(ground_truth.mean), 1, 'first'))];

        all_results.c2l_reg.p_uncorrected = [all_results.c2l_reg.p_uncorrected, mdl_linear.pval];
        all_results.c2l_reg.p_corrected = [all_results.c2l_reg.p_corrected, NaN];
        all_results.c2l_reg.R2ordinary  = [all_results.c2l_reg.R2ordinary, mdl_linear.R2ordinary];
        all_results.c2l_reg.R2adjusted  = [all_results.c2l_reg.R2adjusted, mdl_linear.R2adjusted];
        all_results.c2l_reg.intercept   = [all_results.c2l_reg.intercept, mdl_linear.fit.Coefficients.Estimate(1)];
        all_results.c2l_reg.x1          = [all_results.c2l_reg.x1, mdl_linear.fit.Coefficients.Estimate(2)];
        all_results.c2l_reg.x2          = [all_results.c2l_reg.x2, mdl_linear.fit.Coefficients.Estimate(3)];
        all_results.c2l_reg.opt         = [all_results.c2l_reg.opt, mdl_linear.opt];

        all_results.cluster_based_perm_sig = [all_results.cluster_based_perm_sig, cluster_based_perm_sig];
    end

    if ~exist([path_current.save '\' num2str(idx_subject,formatSpec)], 'dir')
        mkdir([path_current.save '\' num2str(idx_subject,formatSpec)])
    end

end      % for idx_subject = subjects  % iterate over subjects


%% p-value correction

alpha = 0.05;
all_results.c2lcorr.p_corrected = multicmp(all_results.c2lcorr.p_uncorrected,'fdr',alpha);
all_results.c2l_reg.p_corrected = multicmp(all_results.c2l_reg.p_uncorrected,'fdr',alpha);


%% Save variables

save([path_current.save '\Summary_data_allsubjects'], "all_results", "all_trials", '-v7.3');
