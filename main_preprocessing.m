%% main_preprocesing

% Lists the scripts used for preprocessing
% Requires EEGlab and fieldtrip including tesa on the path

%% Initialise EEGlab

eeglab


%% Load the raw data and epoch it

% Loads the raw data
% split up EEG and EMG data
% epoch data

epoch_data


%% preprocess epoched data

preprocessing_pipeline


%% Extract simulation data

simulation_data_extraction


%% Analyse real data

analysis_real_data


%% Calculate Phase accuracy

phase_accuracy

