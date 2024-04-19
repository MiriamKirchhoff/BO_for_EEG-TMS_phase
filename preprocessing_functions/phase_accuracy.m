%% Phase estimation accuracy
%
% Script to assess the phase estimation accuracy of the Lowerlimb data.
% Script adapted from Christoph Zrenner.
%
% version   1.0, 19.04.2024
% author    Miriam Kirchhoff
% project   C2B


clear all

%% Settings

path.data_load = '\rs_EEG';
path.toolboxes = '\Toolboxes';
path.data_save = '\rs_EEG';

settings.participant_names = ...
    {'001','002','003', '004', '005', '007', '008', '009', '011', ...
    '012', '013', '014', '015', '017', '018', '019', '020', '021', ...
    '022', '023', '024', '025', '026', '027', '028', '029', '031', ...
    '032', '034', '035', '036', '038', '039', '040', '041', '042', ...
    '044', '045'};
settings.fs_init = 5000;
settings.fs_ds = 500;
disp(settings)

peak_freq_interval = [9 13];


%% Initialize functions

ang_diff = @(x, y) angle(exp(1i*x)./exp(1i*y));
ang_var = @(x) 1-abs(mean(exp(1i*x)));
design_phastimate_filter = @(ord, freq, fs) ...
    designfilt('bandpassfir', 'FilterOrder', ord, ...
    'CutoffFrequency1', freq-1, 'CutoffFrequency2', freq+1, ...
    'SampleRate', fs, 'DesignMethod', 'window');

eeglab


%% Iterate over participants

for idx_participant = 1:length(settings.participant_names)

    current_participant_name = settings.participant_names{idx_participant};

    %% Load resting state EEG data
    disp(['Start loading data, Subject ' current_participant_name]); tic
    path.rawdata = [path.data_load '\sub-' current_participant_name '_rs_eeg.set'];
    EEG_data = pop_loadset(path.rawdata);

    fprintf('Completed loading data. Elapsed time: %.0f seconds \n', toc)

    EEG_data.label = []
    for i = 1:length(EEG_data.chanlocs)
        EEG_data.label{i} = EEG_data.chanlocs(i).labels;
    end


    %% Filter resting state EEG

    disp('Start C3/C4 Hjorth filter'); tic
    location = 'C3';
    % stimulation on the right hemisphere
    lap_EEG = apply_laplacian_montage(location, EEG_data);
    fprintf('Completed C3/C4 Hjorth filtering. Expired time: %.0f seconds \n', toc)

    lap_EEG = lap_EEG.data;


    %% Calculate and save PSD

    [all_powerspectrum{idx_participant},f_powerspectrum] = ...
        pspectrum(lap_EEG', settings.fs_init, ...
        'FrequencyLimits', [0 45], 'FrequencyResolution', 4);


    %% Downsample

    lap_EEG = resample(lap_EEG, settings.fs_ds, settings.fs_init);


    %% Create different filters

    % Create overlapping epochs from continuous data
    epochs = create_epochs_overlapping(lap_EEG, settings.fs_ds);
    subplot_index = 1;
    ax = subplot(10,5,subplot_index);

    [peak_frequency, peak_SNR{idx_participant}] = ...
        estimate_SNR(epochs, settings.fs_ds, peak_freq_interval, ax);

    if isempty(peak_frequency)
        continue
    end

    % set-up equivalent filter objects for given peak frequency
    filter_objects = {};
    hilbertwin = 128; % this is an appropriate window for alpha at 1000 Hz

    tic
    for ord = [2 3 4 5] % FIR - windowed sinc
        filter_objects = {filter_objects{:} designfilt('bandpassfir', ...
            'FilterOrder', round(ord * (settings.fs_ds/peak_frequency)), ...
            'CutoffFrequency1', peak_frequency-1, ...
            'CutoffFrequency2', peak_frequency+1, ...
            'SampleRate', settings.fs_ds, 'DesignMethod', 'window')};
    end
    toc
    for ord = [3 4 5] % FIR - least squares (equiripple is similar)
        filter_objects = {filter_objects{:} designfilt('bandpassfir', ...
            'FilterOrder', round(ord * (settings.fs_ds/peak_frequency)), ...
            'StopbandFrequency1', peak_frequency-4, ...
            'PassbandFrequency1', peak_frequency-1, ...
            'PassbandFrequency2', peak_frequency+1, ...
            'StopbandFrequency2', peak_frequency+4, ...
            'SampleRate', settings.fs_ds, 'DesignMethod', 'ls')};
    end
    toc
    for ord = [4 8 12] % IIR - butterworth
        filter_objects = {filter_objects{:} designfilt('bandpassiir', ...
            'FilterOrder', ord, 'HalfPowerFrequency1', peak_frequency-1, ...
            'HalfPowerFrequency2', peak_frequency+1, ...
            'SampleRate', settings.fs_ds, 'DesignMethod', 'butter')};
    end
    toc
    for ord = [4 6 8] % IIR - chebychev I
        filter_objects = {filter_objects{:} designfilt('bandpassiir', ...
            'FilterOrder', ord, 'PassbandFrequency1', peak_frequency-1.5, ...
            'PassbandFrequency2', peak_frequency+1.5, ...
            'PassbandRipple', 0.5, 'SampleRate', settings.fs_ds, ...
            'DesignMethod', 'cheby1')};
    end
    toc
    for attenuation = [10 20] % IIR - elliptic
        filter_objects = {filter_objects{:} designfilt('bandpassiir', ...
            'StopbandFrequency1', peak_frequency-2, 'PassbandFrequency1', ...
            peak_frequency-1, 'PassbandFrequency2', peak_frequency+1, ...
            'StopbandFrequency2', peak_frequency+2, ...
            'StopbandAttenuation1', attenuation, ...
            'PassbandRipple', 0.5, 'StopbandAttenuation2', ...
            attenuation, 'SampleRate', settings.fs_ds, ...
            'DesignMethod', 'ellip', 'MatchExactly', 'passband')};
    end
    toc

    [truephase_mean, truephase_variance, trueamp_mean, trueamp_variance] = ...
        phastimate_truephase(epochs, filter_objects);


    %% Settings: Phastimate parameters to be tested for optimization

    filter_order_range = 50:95;     % different filter orders

    bounds_filter_order = [filter_order_range(1) filter_order_range(end)];
    bounds_window = [300 500];      % different window sizes
    bounds_edge = [30 100];         % different edge sizes
    bounds_ar_order = [5 60];       % different orders

    nSamples = ceil(size(epochs, 1)/2);

    % Test whether settings are possible in combination with dataset
    assert(max(bounds_window) <= nSamples)
    assert(min(bounds_window) - 2*max(bounds_edge) > max(bounds_ar_order), ...
        'Parameters can yield empty or too short windows (invalid parameter combinations)')
    assert(min(bounds_window) > 3*max(filter_order_range), ...
        'filtfilt requires data to be at least three times as long as filter order')


    filter_objects_by_order = {}; %the index has to correspond to the order for the genetic algorithm
    for ord = filter_order_range
        filter_objects_by_order{ord} = design_phastimate_filter(ord, peak_frequency, settings.fs_ds);
    end

    % subselect according to true amplitude
    includemask = trueamp_mean >= quantile(trueamp_mean, 0.5);

    [optimal_parameters, ga_output] = phastimate_optimize(epochs(1:ceil(end/2),includemask), truephase_mean(includemask), filter_objects_by_order, bounds_filter_order, bounds_window, bounds_edge, bounds_ar_order, hilbertwin);

    % rerun phastimate with the optimized settings to confirm result
    D = design_phastimate_filter(optimal_parameters.filter_order, peak_frequency, settings.fs_ds);
    [estphase, estamp] = phastimate(epochs(((-optimal_parameters.window_length+1):0)+ceil(end/2),:), D, optimal_parameters.edge, optimal_parameters.ar_order, 128);

    % sanity check if the angular deviation matches the result of the optimization
    phases_error = ang_diff(truephase_mean, estphase);
    assert(abs(optimal_parameters.fval - ang_var(phases_error(includemask))) < 0.01)

    all_optimal_parameters{idx_participant} = optimal_parameters;
    all_phase_error{idx_participant,:} = phases_error;
    all_meanphase_error(idx_participant) = circ_mean(phases_error);

end


% Save data

disp('Start saving data'); tic
save([path.data_save 'EEG_phase_accuracy.mat'], 'all_phase_error', 'all_meanphase_error', 'peak_SNR', 'all_optimal_parameters', 'all_powerspectrum', 'f_powerspectrum', '-v7.3')
fprintf('Completed saving data. Elapsed time: %.0f seconds \n', toc)


%% Calculate SNR

for idx_participant = 1:length(settings.participant_names)

    current_participant_name = settings.participant_names{idx_participant};

    %% Load resting state EEG data
    disp(['Start loading data, Subject ' current_participant_name]); tic
    path.rawdata = [path.data_load '\sub-' current_participant_name '_rs_eeg.set'];
    EEG_data = pop_loadset(path.rawdata);
    %clear data
    fprintf('Completed loading data. Elapsed time: %.0f seconds \n', toc)

    EEG_data.label = [];
    for i = 1:length(EEG_data.chanlocs)
        EEG_data.label{i} = EEG_data.chanlocs(i).labels;
    end


    %% Filter resting state EEG

    disp('Start C3/C4 Hjorth filter'); tic
    location = 'C3';
    % stimulation on the right hemisphere
    lap_EEG = apply_laplacian_montage(location, EEG_data);
    fprintf('Completed C3/C4 Hjorth filtering. Expired time: %.0f seconds \n', toc)

    lap_EEG = lap_EEG.data;


    SNR(idx_participant) = snr(lap_EEG, 5000, 5);

    lap_EEG_bp = bandpass(lap_EEG,[9 13],5000);

    SNR_bp(idx_participant) = snr(lap_EEG_bp, 5000, 5);

end

%% Get means and std

sig_participants = ismember(settings.participant_names, {'001' '007' '019' '022' '032'});
mean_SNR = mean(SNR)
std_SNR = std(SNR)

mean_SNR_sig = mean(SNR(sig_participants))
std_SNR_sig = std(SNR(sig_participants))

% Get all phase errors
errors_all = [];
for idx_participant = 1:length(settings.participant_names)
    errors_all = [errors_all all_phase_error{idx_participant,:}];
end

mean_errors = rad2deg(mean(errors_all))
std_errors = rad2deg(std(errors_all))

mean_errors_sig = rad2deg(mean(errors_all(sig_participants)))
std_errors_sig = rad2deg(std(errors_all(sig_participants)))

%% Save data

save([path.save 'phase_acc_summary']);
