function  results = mleBIC(pos,models,R,dt)
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

% BIC model Selection
numData = (size(pos,1)-1)*size(pos,2);
nParams = zeros(length(models),1);
for i = 1:length(models)
    nParams(i) = length(mu{i});
end
BIC = maxLogL - nParams/2*log(numData);   
modelProb = exp(BIC - max(BIC));
modelProb = modelProb'/sum(modelProb);

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


