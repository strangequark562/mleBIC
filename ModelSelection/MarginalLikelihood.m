function modelProb = MarginalLikelihood(maxLogL,mu,I,strength)
%-------------------------------------------------------------------------- 
% Summary: This function calculates the marginal likelihood for each model
% using the Laplace approximation.
% 
% Input:
%       maxLogL = vector of maximum log-likelihood values for each model
%       mu = cell of vectors of MLE estimates for each model
%       I = cell of observed Fisher information matrices for each model
%       strength = tuning parameter for priors
%
% Output:
%       modelProb = probability that trajectory generated from each model
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

numModels = length(maxLogL);

% calculate logEvidence
logPrior = zeros(numModels,1);
curvature = zeros(numModels,1);
for j = 1:numModels
    logPrior(j) = LogPrior(mu{j},mu{j},I{j}*strength);
    curvature(j) = 0.5*length(mu{j})*log(2*pi) + 0.5*LogDeterminant(I{j});
end
logEvidence = maxLogL + logPrior + curvature;

% normalize evidence
modelProb = exp(logEvidence - max(logEvidence));
modelProb = modelProb'/sum(modelProb);
