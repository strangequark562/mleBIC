function logL = Likelihood_DrivenDiffusion(N,deltax,b,R,dt)
%-------------------------------------------------------------------------- 
% Summary: Likelihood model for trajectories undergoing driven diffusion
% 
% Input:
%       N = number of displacements
%       deltax = vector of x-displacements
%       deltay = vector of y-displacements
%       D = diffusion coefficient 
%       v = vector of drift velocities, v = [vx vy];
%       sigma = static localization noise
%       R = motion blur coefficient
%       dt = frame duration (s)
%
% Output:
%       logL = log-likelihood value for a given set of model parameters
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

D = b(1);
v = b(2:end-1);
sigma = b(end);

% covariance matrix
msd = 2*D*dt;
staticnoise = 2*sigma^2;
motionblur = -2*R*msd;
C = toeplitz([msd+staticnoise+motionblur -(staticnoise+motionblur)/2 zeros(1,N-2)]);
 
% log-likelihood
mu = v*dt;
logL = LogMultiVariateNormal(N,deltax,C,mu);

end
