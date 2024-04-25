function [sample, subset_mean, subset_sd, subset_ci, N] = simulate_data_slidewindow(X, dataset, wintype, winsize, limits, robust_estimate, surrogate)
% function [sample, subset_mean, subset_sd, subset_ci] = simulate_data_slidewindow(X, dataset, wintype, winsize, limits)
%
% Draws a sample at location X based on a data set.
% Uses a sliding window around point X to define a normal distribution with
% mean and sd from which to draw the sample.
%
% INPUTS
    % X:        X-value at which is being sampled, double
    % dataset:  Dataset from which samples are being drawn, samplesize*2 matrix with
    %           x-values in first column and y-values in second column
    % wintype:  define whether the window is limited by number of samples
    %           ('sample_absolute'), percentage of total samples in decimals
    %           ('sample_relative'), or total windowsize ('winsize'), string
    % winsize:  define window size based on windowtype
    % limits:   define min and max period of the data as vector [min max]
    % surrogate:number of points to be calculated, CI will be bonferroni corrected 
    % robust_estimate: If true, replace mean and std by median and mad,
    %           default is false, bool.
%
% OUTPUTS
    % sample:   Y-value sampled from dataset corresponding to X, double
    % subset_mean: mean of the normal distribution from which was sampled,
    %           double
    % subset_sd:SD of the normal distribution from which was sampled,
    %           double
    % subset_ci: 95% CI of the normal distribution from which was sampled,
    %           double
%
% version   1.1, 19.04.2023
% author    Miriam Kirchhoff
% project   C2B

if nargin<6
    robust_estimate = false;
end

if nargin<7
    surrogate = 1;
end

% Account for periodicity of the data: copy dataset to p-1 and p+1
period = max(limits) - min(limits);
dataset_lower = dataset;
dataset_lower(:,1) = dataset(:,1) - period;
dataset_higher = dataset;
dataset_higher(:,1) = dataset(:,1) + period;
dataset_extended = [dataset;dataset_lower;dataset_higher];

% Define relevant subset of data from dataset for fitting gaussian distribution
switch wintype
    case {'sample_absolute', 'sample_relative'}
        if strcmp(wintype, 'sample_relative') % calculate sample size based on selected method
            samplesize = round(winsize*size(dataset,1));
        elseif strcmp(wintype, 'sample_absolute')
            samplesize = winsize;
        end % if strcmp(wintype, 'sample_relative') % calculate sample size based on selected method

        % sort data by distance to X
        distances = abs(dataset_extended(:,1) - X);
        [~, sort_idx] = sort(distances);

        % select #samplesize closest samples to be in subset
        relevant_locations = sort_idx(1:samplesize);

    case 'winsize'
        % borders defined as X+-0.5*winsize
        min_X = X-0.5*winsize;
        max_X = X+0.5*winsize;

        % relevant data has X values between these borders
        relevant_locations = dataset_extended(:,1) > min_X & dataset_extended(:,1) < max_X;
end % switch wintype

% select rows that meet criterion
subset = dataset_extended(relevant_locations,:);

% use the subset of data to fit gaussian
if robust_estimate
    subset_mean = median(subset(:,2));
    subset_sd = mad(subset(:,2));
else
    subset_mean = mean(subset(:,2));
    subset_sd = std(subset(:,2));
end

N = length(subset(:,2));

alpha = 0.05/surrogate;

SEM = subset_sd/sqrt(length(subset(:,2)));               % Standard Error
ts = tinv([0.025  0.975],length(subset(:,2))-1);      % T-Score
subset_ci = subset_mean + ts'*SEM;                      % Confidence Intervals   

% draw a random sample from gaussian
sample = normrnd(subset_mean,subset_sd);

end