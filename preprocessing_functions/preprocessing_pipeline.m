%% preprocessing_pipeline
%
% Preprocess epoched data after epoch_data for the BO simulation
%
% version   1.0, 19.04.2024
% author    Miriam Kirchhoff
% project   C2B

clear
close all
clc

%% Settings

subjects = 1:38;
formatSpec = '%03.0f';

% settings
path.load = '\data_epoched\';
path.save = '\data_preprocessed\';
mkdir(data.save)

% settings for montage
hjorth.channel = 'C3';
switch hjorth.channel  % get surrounding variables for hjorth filter
    case 'C3'
        hjorth.electrodes = {'FCC5h', 'FCC3h', 'CCP5h', 'CCP3h'};
end

% settings for downsampling
ds_fs = 1000;       % downsampling freqency in Hz

%% Iterate over subjects

for idx_subject = subjects  % iterate over subjects

    clearvars -except idx_subject path hjorth formatSpec


    %% Load epoched data

    % start timer
    tic

    % start notification
    fprintf('\nStarting participant %03.0f. \n', idx_subject)
    disp('Start loading...')

    % load EEG data
    load([path.load '\' num2str(idx_subject,formatSpec) '\EEG_P' num2str(idx_subject,formatSpec)]);
    EEG = EEG_epoched;
    clear EEG_epoched

    % timing notification
    fprintf('loading completed. Expired time: %.0f seconds \n', toc)


    %% Find Hjorth electrodes

    temp = struct2cell(EEG.chanlocs);
    hjorth.channel_idx = find(strcmpi(hjorth.channel, temp(1,:)));
    for i = 1:length(hjorth.electrodes)
        hjorth.electrode_idx(i) = find(strcmpi(hjorth.electrodes{i}, temp(1,:)));
    end

    clear i


    %% Baseline correction
    % Subtract baseline from each channel location and epoch

    % settings
    % define prestim period
    t.baseline = [-600 -5];         % ms time wrt TMS

    % start timer
    tic

    % start notification
    disp('Start baseline correction...')

    % remove baseline
    EEG = mk_rmbase(EEG, t.baseline);

    % timing notification
    fprintf('baseline correction completed. Expired time: %.0f seconds \n', toc)


    %% Detrend

    % settings
    % define prestim period
    t.detrend = [min(EEG.times), min(-5, max(EEG.times))];  % ms time wrt TMS

    % start timer
    tic

    % start notification
    disp('Start detrending...')

    % detrend data linearly
    EEG = tesa_detrend(EEG, 'linear', t.detrend);

    % timing notification
    fprintf('detrending completed. Expired time: %.0f seconds \n', toc)


    %% Find bad channels and exclude trials where the threshold was crossed for any channel in the laplacian

    % settings
    % define time relevant for exclusion
    t.channel_threshold = [-600 -5];
    rejection.threshold_level = 150;

    % start timer
    tic

    % start notification
    disp('Start channel rejection by threshold...')

    % find relevant channels per trial
    threshold_marker = range(EEG.data ...
        (:, EEG.times >= t.channel_threshold(1) & ...
        EEG.times <= t.channel_threshold(2),:),2)>rejection.threshold_level;

    % find trials that cross the threshold in any relevant channel for the
    % montage
    temp = squeeze(threshold_marker([hjorth.electrode_idx hjorth.channel_idx],:,:));

    % note trials where the montage is not reliable
    rejection.threshold = sum(temp) > 1;
    clear temp
    fprintf('%.0f trials rejected by threshold. \n', sum(rejection.threshold))

    % reject selected trials
    EEG = pop_select(EEG, 'trial', find(~rejection.threshold));

    % timing notification
    fprintf('threshold rejection completed. Expired time: %.0f seconds \n', toc)


    %% Downsample

    % start timer
    tic

    % start notification
    disp('Start downsampling...')

    % downsample data
    EEG = pop_resample(EEG, ds_fs);

    % timing notification
    fprintf('downsampling completed. Expired time: %.0f seconds \n', toc)


    %% BP filter

    % settings
    bp.frequencies = [9 13];

    % start timer
    tic

    % start notification
    disp('Start bp filtering...')

    % bp filter
    EEG = pop_eegfiltnew(EEG, 'locutoff', bp.frequencies(1),'hicutoff', bp.frequencies(2));

    % timing notification
    fprintf('BP filtering completed. Expired time: %.0f seconds \n', toc)


    %% Laplacian montage

    % start timer
    tic

    % start notification
    disp('Start laplacian...')

    % define laplacian montage as filter
    lap_filter = zeros(size(EEG.data,1),1);
    lap_filter(hjorth.channel_idx) = 1;
    lap_filter(hjorth.electrode_idx) = -1/4;

    EEG_lap = EEG;

    % apply laplacian montage
    temp = (reshape(EEG_lap.data, 126, [])'*lap_filter)';
    temp = reshape(temp, length(EEG_lap.times), []);
    EEG_lap.data = temp;
    EEG_lap.chanlocs = 'C3-centered laplacian';
    EEG_lap.nbchan = 1;
    EEG_lap.soundLeadfield = [];

    % timing notification
    fprintf('Laplacian montage completed. Expired time: %.0f seconds \n', toc)


    %% Save final data

    % start timer
    tic

    % start notification
    disp('Start saving data...')

    % create folder if not existent
    if ~exist([path.save '\' num2str(idx_subject,formatSpec)], 'dir')
        mkdir([path.save '\' num2str(idx_subject,formatSpec)])
    end

    % save data
    save([path.save '\' num2str(idx_subject,formatSpec) '\EEG_P' num2str(idx_subject,formatSpec)], "EEG", "t", "path", "bp", "rejection", '-v7.3');
    save([path.save '\' num2str(idx_subject,formatSpec) '\EEG_lap_P' num2str(idx_subject,formatSpec)], "EEG_lap", "t", "path", "bp", "rejection", '-v7.3');


    % timing notification
    fprintf('Saving completed. Expired time: %.0f seconds \n', toc)

end     % for idx_subject = subjects  % iterate over subjects






