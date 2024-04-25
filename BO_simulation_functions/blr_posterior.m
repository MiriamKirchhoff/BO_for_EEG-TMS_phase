function [m_N, S_N, S_N_inv] = blr_posterior(phi, t, alpha, beta)
% function [m_N, S_N, S_N_inv] = blr_posterior(phi, t, alpha, beta)
%
% Computes mean and covariance matrix of the posterior distribution.
%
% INPUTS
    % phi:  double n*(d+1)  design matrix phi(X), size is
    %                       num_samples*regression dimensions + intercept
    % t:    double n*1      y-values of past measurements
    % alpha:scalar          weight variance of regression
    % beta: scalar          noise variance of regression
% OUTPUTS    
    % m_N   double (d+1)*1  mean vector of posterior, dimensions +
    %                       intercept
    % S_N   double (d+1)*(d+1) cov matrix of posterior
    % S_N_inv double (d+1)*(d+1) inverse cov matrix of posterior
%
% version   1.0, 04.08.2023
% author    Miriam Kirchhoff, adapted from 
% project   C2B

% calculate inverse
S_N_inv = alpha * eye(size(phi, 2)) + beta * phi' * phi;
S_N = inv(S_N_inv);

% calculate mean
m_N = beta * S_N * phi' * t;

end
