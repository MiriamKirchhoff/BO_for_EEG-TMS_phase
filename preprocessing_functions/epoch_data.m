%% epoch_data
%
% Script to epoch and merge the raw dataset. Data is being saved
% separately as EMG and EEG data with an epoch around the stimulation onset
% for data size reduction proior to the following data preprocessing. 
%
% version   1.0, 19.04.2024
% author    Miriam Kirchhoff
% project   C2B

clear all

% subject indices
loaddata.subjects = 1:38;

% define epoch windows
epoch_window_EEG = [-.704 -0.005];      % ms
epoch_window_EMG = [-.1 .1];            % ms

% names of the EMG channels
mep_amp.channels = {'APBr', 'FDIr'};    % names of MEP channels
formatSpec = '%03.0f';

% paths
path.load = '\sub-';
path.save = '\data_epoched';

% create save location
mkdir(path.save)


for idx_subject = loaddata.subjects   % iterate over subjects

    % Participant start information
    fprintf('\nStarting participant %03.0f\n', idx_subject)


    %% Load datafile

    path.current = [path.load_current num2str(idx_subject, formatSpec) ...
        '\eeg\sub-' num2str(idx_subject, formatSpec) '_eeg.set'];
    EEG_raw = pop_loadset(path.current);
    fprintf('loading completed. Expired time: %.0f seconds \n', toc)


    %% Split dataset into EMG and EEG

    % section timer
    tic

    % find EMG channel locations
    temp = struct2cell(EEG_raw.chanlocs);
    for i = 1:length(mep_amp.channels)
        mep_amp.channel_idx(i) = find(strcmpi(mep_amp.channels{i}, temp(1,:)));
    end % for i = 1:length(mep_amp.channels)

    % split
    EEG = pop_select(EEG_raw, 'rmchannel', mep_amp.channel_idx);
    EMG = pop_select(EEG_raw, 'channel', mep_amp.channel_idx);

    fprintf('splitting dataset completed. Expired time: %.0f seconds \n', toc)

    clear EEG_raw


    %% Epoch data

    % section timer
    tic

    % Epoch data
    [EEG_epoched, idx] = pop_epoch(EEG, 'A - Stimulation', epoch_window_EEG);
    [EMG_epoched, idx] = pop_epoch(EMG, 'A - Stimulation', epoch_window_EMG);
    clear idx

    fprintf('epoching completed. Expired time: %.0f seconds \n', toc)

    clear EEG EMG


    %% Save data of one participant

    % section timer
    tic

    % create folder if not existent
    if ~exist([path.save '\' num2str(idx_subject,formatSpec)], 'dir')
        mkdir([path.save '\' num2str(idx_subject,formatSpec)])
    end

    % save data
    save([path.save '\' num2str(idx_subject,formatSpec) '\EEG_P' num2str(idx_subject,formatSpec)], 'EEG_epoched', '-v7.3');
    save([path.save '\' num2str(idx_subject,formatSpec) '\EMG_P' num2str(idx_subject,formatSpec)], 'EMG_epoched', '-v7.3');

    fprintf('saving completed. Expired time: %.0f seconds \n', toc)

    clear EEG_epoched EMG_epoched


end % for idx_subject = 1:length(loaddata.subjects) % iterate over subjects

