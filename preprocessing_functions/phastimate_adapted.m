function [phase, amplitude] = phastimate_adapted(data, edge, ord, hilbertwindow, varargin)
%PHASTIMATE Pre-stimulus phase estimation
%   phase = PHASTIMATE(data, D, edge, ord, hilbertwindow, [offset_correction], [iterations], [armethod]) estimates the phase of the matrix data
%
%   data is a time x channel matrix
%   D is the digitalFilter object
%   edge is the filter edge size (removed)
%   ord is the order of the autoregressive yule-walker model
%   hilbertwindow is the length of the window used for the hilberttransform
%   iterations is the number of samples forward predicted (default is edge + hilbertwindow/2)
%   armethod is the method used to estimate autoregressive parametres (default is @aryule)
%
%   Example:
%   Fs = 200;
%   data = [sin([1:Fs]/Fs*5*2*pi)' sin([1:Fs]/Fs*6*2*pi)'] + randn(Fs,2);
%   D = designfilt('bandpassfir', 'FilterOrder', round(Fs/5), 'CutoffFrequency1', 4, 'CutoffFrequency2', 7, 'SampleRate', Fs);
%
%   phase = phastimate(data, D, 25, 20, 64)
%
%   Adapted from Christoph Zrenner by Miriam Kirchhoff

if nargin < 5
    offset_correction = 0;
else
    offset_correction = varargin{1};
end
if nargin < 6
    iterations = edge + ceil(hilbertwindow/2);
else
    iterations = varargin{2};
    assert(iterations > edge, 'iterations must be larger than the number of edge samples')
end
if nargin < 7
    armethod = @aryule; %could be aryule, arburg
else
    armethod = varargin{3};
end


% demean the data
data = detrend(data,'constant');

% cut the edge
data_filtered_withoutedge = data(edge+1:end-edge,:);

% determine AR parameters
[a, e, rc] = armethod(data_filtered_withoutedge, ord);
coefficients = -1 * flip(a(:, 2:end)');

% prepare matrix with the aditional time points for the forward prediction
data_filtered_withoutedge_predicted = [data_filtered_withoutedge; ones(iterations, size(data,2))];
% run the forward prediction
for i = iterations:-1:1
    data_filtered_withoutedge_predicted(end-i+1,:) = ...
        sum(coefficients .* data_filtered_withoutedge_predicted((end-i-ord+1):(end-i),:));
end

data_filtered_withoutedge_predicted_hilbertwindow = data_filtered_withoutedge_predicted(end-hilbertwindow+1:end,:);

% analytic signal and angle
data_filtered_withoutedge_predicted_hilbertwindow_analytic = hilbert(data_filtered_withoutedge_predicted_hilbertwindow);

phase = angle(data_filtered_withoutedge_predicted_hilbertwindow_analytic(end-iterations+edge+offset_correction,:));
amplitude = mean(abs(data_filtered_withoutedge_predicted_hilbertwindow_analytic));

end