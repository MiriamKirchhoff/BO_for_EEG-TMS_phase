function [xnext, KG] = blr_nextpoint_kg(X, t, beta, m_N, S_N, T)
% function [xnext, KG] = blr_nextpoint_kg(X, t, beta, m_N, S_N, T)
%
% Calculates the next point as well as the acquisition function for a
% circular to linear regression function using knowledge gradient. 
%
% INPUTS
    % X:    double n*1      x-values of past measurements
    % t:    double n*1      y-values of past measurements
    % beta: scalar          noise variance of regression
    % m_N   double (d+1)*1  mean vector of regression, dimensions +
    %                       intercept
    % S_N   double (d+1)*(d+1) cov matrix of regression
    % T     double N*1      possible measurement locations
% OUTPUTS
    % xnext:scalar          optimal next measurement location
    % KG:   double N*1      acquisition function for possible measurement
    %                       locations
%
% version   1.0, 04.08.2023
% author    Miriam Kirchhoff, adapted from J Sarvas, A Tervo (2020)
% project   C2B

uselog = 1;
n = length(X);
N = length(T);
% N=length(T);

% Compute covariances between measurement points
phi_X = [ones(size(X)), sin(X), cos(X)];
G = 1/beta + phi_X * S_N * phi_X';
mu_X = blr_posterior_predictive(phi_X, m_N, S_N, beta);
w = (G + eye(n)/beta)\(t - mu_X); % same as w=invGL*(t-mu0) but faster and more accurate

% Compute covariances between possible measurement points
phi_T = [ones(size(T)), sin(T), cos(T)];
mu_T = blr_posterior_predictive(phi_T, m_N, S_N, beta);
H1 = 1/beta + phi_T * S_N * phi_T';

% Compute covariances between possible measurement points and measurement
% points
H2 = 1/beta + phi_T * S_N * phi_X';

% Compute the knowledge gradient
S = H1 - (H2*((G + eye(n)/beta)\H2'));

for j=1:N
    a = mu_T + H2 * w;
    b = S(:,j)./sqrt(S(j,j) + beta);
    [a,b]=tms_prep(a,b);
    [c,J]=tms_breakpoints(a,b);
    aa=a(J);
    bb=b(J);

    if uselog % Compute KG on log scale (better numerical accuracy)
        s = -abs(c);
        logsumx2 = zeros(size(s));
        inds = find(s<-10); % numerically unstable region
        if numel(inds)>0
            logsumx2(inds) = -0.5*log(2*pi)-s(inds).^2/2-log(s(inds).^2+1);
        end
        inds = find(s>=-10); % numerically stable region
        if numel(inds)>0
            logsumx2(inds) = log(normpdf(s(inds))+s(inds).*normcdf(s(inds)));
        end
        logsumx = log(diff(bb))'+logsumx2;
        En = max(logsumx)+log(sum(exp(logsumx-max(logsumx))));
    else % KG in linear scale (numerically less stable)
        s=0;
        for i=1:length(c)
            u=c(i);
            s=s+(bb(i+1)-bb(i))*(normpdf(u)+u*normcdf(u));
        end
        En=s+aa(end);
    end
    if En
        KG(j)=En;
    else
       KG(j)=NaN;
    end
end
maxKG=max(KG);
k = find(KG==maxKG);
if numel(k)>1 % If several points with max value, pick one of them randomly
    k = k(randi(numel(k)));
end
xnext=k; % Return index of KG max


% SUBROUTINES
% As in Frazier, P., Powell, W., Dayanik, S., 2009. The knowledge-gradient policy 
% for correlated normal beliefs. Inf. J. Comput. 21, 599–613. 
% https://doi.org/10.1287/ ijoc.1080.0314.
 
function [u,v]=tms_prep(a,b)
[b,ind]=sort(b);
a=a(ind);
d=diff(b);
n=length(a);

u=a(1);
v=b(1);
for i=2:n
    if b(i)-v(end)>1e-10
        u=[u;a(i)];
        v=[v;b(i)];
    else
        u(end)=max(u(end),a(i));
    end
end

function [c,J]=tms_breakpoints(a,b)
m=length(a);
c=[];
J=1;
i=1;
while i<m
    ind=i+1:m;
    [s,k]=min((a(i)-a(ind))./(b(ind)-b(i)));
    c=[c,s];
    J=[J,i+k];
    i=i+k;
end



