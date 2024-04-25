function [m_N, S_N, beta] = blr_regression(X, t)
% Uses a linear basis function bayesian regression to fit a sine function
% to data X, t(X)
% m_n(1) + m_n(2)*sin(X) + m_n(3)*cos(X)
% y(X) = N(m_n, S_N)
%
% INPUTS
    % X:    double n*1      x-values of past measurements
    % t:    double n*1      y-values of past measurements
%
% OUTPUTS    
    % m_N   double (d+1)*1  mean vector of posterior, dimensions +
    %                       intercept
    % S_N   double (d+1)*(d+1) cov matrix of posterior
    % beta: scalar          noise variance of regression
%
% version   1.0, 04.08.2023
% author    Miriam Kirchhoff, adapted from J Sarvas, A Tervo (2020)
% project   C2B


%% Settings

% set initial values for alpha and beta
alpha_0 = 1e-10;
beta_0  = 1e-10;

% max number of iterations
max_iter = 600;

% convergence criterion
rtol = 1e-10;


%% Bayesian linear regression

% calculate design matrix
phi = [ones(size(X)), sin(X), cos(X)];

% initialise
eigenvalues_0 = eig(phi'*phi);
beta = beta_0;
alpha = alpha_0;

% define number of samples and number of weights
[N, ~] = size(phi);

for i = 1:max_iter

    % save previous alpha for convergence criterion
    alpha_prev = alpha;

    % scale eigenvalues
    eigenvalues = (eigenvalues_0 * beta) + 1e-5;

    % compute mean and covariance matrix of posterior distribution
    [m_N, S_N, S_N_inv] = blr_posterior(phi, t, alpha, beta);

    % update alpha and beta
    gamma = sum(eigenvalues./(eigenvalues + alpha));
    alpha = gamma / sum(m_N.^2);
    beta_inv = 1 / (N - gamma) * sum((t - phi * m_N) .^ 2);
    beta = 1 / beta_inv;

    % check for convergence
    if abs(alpha_prev-alpha) <= rtol
        %fprintf('Convergence after %d iterations.\n', i + 1)
        break
    end

end
