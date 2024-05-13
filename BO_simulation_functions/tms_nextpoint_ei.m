function [xnext, EI]=tms_nextpoint_ei(postmean, postvar, postmean2)

% Function finds the next measurement point by the EI algorithm


best_mu = max(postmean2);
postSD = sqrt(postvar);
z = (postmean - best_mu) ./ postSD;

EI = (postmean - best_mu) .* normcdf(z) + postSD .* normpdf(z);

EI(postSD == 0) = 0;

xnext = find(EI==max(EI));

% if more than one optimum: pick random one
ri = randi(length(xnext), 1);
% plot(postmean)
% hold on
% plot(postmean+postvar)
% plot(postmean-postvar)
% plot(EI*100)
xnext = xnext(ri);











