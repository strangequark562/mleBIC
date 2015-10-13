function logL = Likelihood_Immobile(N,deltax,sigma)
%-------------------------------------------------------------------------- 
% Summary: Likelihood model for immobile trajectories
% 
% Input:
%       N = number of displacements
%       deltax = vector of x-displacements
%       deltay = vector of y-displacements
%       sigma = static localization noise
%
% Output:
%       logL = log-likelihood value for a given sigma
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

% generate covariance matrix
staticnoise = 2*sigma^2;
C = toeplitz([staticnoise -staticnoise/2 zeros(1,N-2)]);

% log-likelihood
logL = LogMultiVariateNormal(N,deltax,C);

end