function logL = Likelihood_confinedDiffusion(N,deltax,D,L,sigma,R,dt)
%-------------------------------------------------------------------------- 
% Summary: Likelihood model for trajectories undergoing confined diffusion
% 
% Input:
%       N = number of displacements
%       deltax = vector of x-displacements
%       deltay = vector of y-displacements
%       D = diffusion coefficient 
%       L = confinement size
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


% diagonals of covariance matrix
 k = 1:2:5000;
msd = L^2/6 - sum(16*L^2/pi^4./k.^4.*exp(-k.^2*pi^2/L^2*D*dt));

% 1st off-diagonal of covariance matrix
corr_1 = -L^2/12+16*L^2/pi^4*sum(1./k.^4.*exp(-(k*pi/L).^2*D*dt))-8*L^2/pi^4*sum(1./k.^4.*exp(-2*(k*pi/L).^2*D*dt));

% higher order off-diagonals of covariance matrix
k = 1:2:5001;
corr_n = zeros(1,N-1);
for n = 1:N-1
	 corr_n(n) = -8*L^2/pi^4*sum(1./k.^4.*(exp(-(k*pi/L).^2*D*n*dt) - 2*exp(-(k*pi/L).^2*D*(1+n)*dt) + exp(-(k*pi/L).^2*D*(2+n)*dt)));
end

% generate noiseless covariance matrix
C_confined = toeplitz([msd corr_1 corr_n]);

% motion blur corrections
vacf = CalculateVACF(C_confined);
motionblur1 = -R*(2*vacf(1) - 2*vacf(2));
motionblur2 = -R*(2*vacf(2:end-1) - vacf(3:end) - vacf(1:end-2) );
C_dynamic = toeplitz([motionblur1 motionblur2]);

% corrections for static localization noise
staticnoise = 2*sigma^2;
C_static = toeplitz([staticnoise -staticnoise/2 zeros(1,N-2)]);

% total covariance matrix
C = C_confined(1:end-1,1:end-1) + C_dynamic + C_static;

% calculate log-likelihood
logL = LogMultiVariateNormal(N,deltax,C);

