function [y, y_var, y_cov] = blr_posterior_predictive(phi, m_N, S_N, beta)
% function [y, y_var, y_cov] = blr_posterior_predictive(phi, m_N, S_N, beta)
%
% Computes mean and variances of the posterior predictive distribution.
%
% INPUTS
    % phi:  double n*(d+1)  design matrix phi(X)
    % m_N   double (d+1)*1  mean vector of posterior, dimensions +
    %                       intercept
    % S_N   double (d+1)*(d+1) cov matrix of posterior
    % beta: scalar          noise variance of regression
% OUTPUTS    
    % y     double n*1      predicted mean at location X
    % y_var double n*1      predicted variance at location X
    % y_cov double n*n      predicted covariance of measurements X * X
%
% version   1.0, 04.08.2023
% author    Miriam Kirchhoff, adapted from J Sarvas, A Tervo (2020)
% project   C2B

y = phi * m_N;
% Only compute variances (diagonal elements of covariance matrix)
y_var = 1 / beta + sum((phi * S_N) .* phi, 2);

% calculate covariance
y_cov = 1/beta + phi * S_N * phi';

end
