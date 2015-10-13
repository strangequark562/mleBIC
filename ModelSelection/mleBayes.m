function  results = mleBayes(pos,models,R,dt,strength)
%-------------------------------------------------------------------------- 
% Summary: mleBayes performs MLE estimates and model selection for a 2D
% trajectory.
% 
% Input:
%       x = vector of x-positions
%       y = vector of y-positions
%       models = diffusion models to analyze
%       R = motion blur coefficient
%       dt = frame duration (s)
%       strength = parameter to tune prior strength
%
% Output:
%       results = structure containing: 
%         results.model = optimal model
%         results.modelProb = optimal model probability
%         results.mu = optimal model parameters
%         results.I = optimal model observed Fisher information
%         results.estimates = same values for other models

% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 


% initial estimates and bounds for model parameters
params = InitialParameters(models,pos,R,dt);

% numerically find maximum log-likelihood estimates
[mu,I,maxLogL] = MaximumLogLikelihood(pos,params);

% calculate marginal likelihood using Laplace Approximation
modelProb = MarginalLikelihood(maxLogL,mu,I,strength);


% find optimal model
[MAX,index] = max(modelProb);
results = struct;
results.model = models{index};
results.modelProb = MAX;
results.mu = mu{index};
results.I = I{index};
results.estimates.I = I;
results.estimates.maxLogL = maxLogL;
results.estimates.modelProb = modelProb;
results.estimates.mu = mu;


