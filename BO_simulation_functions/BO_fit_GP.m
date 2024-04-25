function [y, y_var, GP_fit] = BO_fit_GP(X, t, T, GP_fit_before)
% function [y, y_var, GP_fit] = BO_fit_GP(X, t, T, GP_fit_before)
% fit a gaussian process to the data based on a prior process
%
% INPUTS:
    % X:    double n*1      x-value indices of T of past measurements
    % t:    double n*1      y-values of past measurements
    % T:    double N*1      x-values of potential measurement locations
    % GP_fit_before struct  struct containing the prior model containing
    %                       fields:
        % postmean      double N*1      mean of prior model fit on T
        % postvar       double N*1      variance of prior model fit on T
        % postmeanX     double N*1      mean of prior model fit on T(X)
        % kernel        string          definition of the kernel type used
%
% OUTPUTS:
    % y:    double N*1      mean of posterior model fit on T
    % y_var:double N*1      variance of the posterior model fit on T
    % GP_fit: struct        struct containing posterior model containing
    %                       fields:
        % postmean      double N*1      mean of posterior model fit on T
        % postvar       double N*1      variance of posterior model fit on T
        % postmeanX     double N*1      mean of posterior model fit on T(X)
        % kernel        string          definition of the kernel type used
        % mu0           scalar          column vector of the prior mean of 
        %                               [f(X(:,1)),...,f(X(:,n))]'
        % a0            scalar          noise variance for kernel
        %                               computation
        % a1            scalar          variance for kernel computation
        % lam2          double n*1      column vectors of the noise 
        %                               variances at points in X
        % muz3          scalar          (N,1) column vector of the prior 
        %                               mean of [f(T(:,1)),...,f(T(:,N))]'
%
% version   1.0, 04.08.2023
% author    Miriam Kirchhoff, adapted from J Sarvas, A Tervo (2020)
% project   C2B


% initialization
postmean = GP_fit_before.postmean;             % posterior mean array,
postvar = GP_fit_before.postvar;               % posterior variance array,
postmeanX = GP_fit_before.postmeanX;           % and posterior mean in measurement points array
k = length(postmean);  % define over size      % Zero iteration counter
kernel = GP_fit_before.kernel;

[mu0,muz2,muz3,a0,a1,lam2] = set_GP_parameters(k,t,T,T,postmean,postvar,postmeanX,X);

% Compute posterior mean in computation grid
[postmean,postvar]=tms_postmean(T(X),t',T,mu0,muz2,a0,a1,lam2,kernel,1);

% Compute posterior mean in measurement grid (needed for the parameter
% setting on the next iteration)
postmeanX = tms_postmean(T(X),t',T,mu0,muz3,a0,a1,lam2,kernel,1);

GP_fit.mu0      = mu0;
GP_fit.a0       = a0;
GP_fit.a1       = a1;
GP_fit.lam2     = lam2;
GP_fit.muz3     = muz3;
GP_fit.kernel   = kernel;
GP_fit.postmeanX= postmeanX;
GP_fit.postmean = postmean;
GP_fit.postvar  = postvar;

y       = postmean;
y_var   = postvar;

end

