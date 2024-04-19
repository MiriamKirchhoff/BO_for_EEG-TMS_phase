function EEG_rm = mk_rmbase(EEG, timerange)
% function EEG_rm = mk_rmbase(EEG)
% Removes the baseline of epoched data foch each channel for each trial and
% does not take hours like the EEGlab function
%
% Inputs:
%   EEG        - Input dataset
%   timerange  - [min_ms max_ms] Baseline latency range in milliseconds.
%
% Outputs:
%   EEG_rm      Output dataset with baseline removed

% find time indices that are relevant
t_idx = EEG.times >= timerange(1) & EEG.times <= timerange(2);

% calculate averages within this range for each epoch and channel across
% time
data = EEG.data;
baselines = mean(data(:,t_idx,:),2);

% subtract baselines from data
data_rm = data - baselines;

% insert data into return struct
EEG_rm = EEG;
EEG_rm.data = data_rm;

end