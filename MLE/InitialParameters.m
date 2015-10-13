function params = InitialParameters(models,pos,R,dt)
%-------------------------------------------------------------------------- 
% Summary: Estimate initial parameters based on trajectory properties
% 
% Input:
%       models = diffusion models to analyze
%       x = vector of x-positions
%       y = vector of y-positions
%       R = motion blur coefficient
%       dt = frame duration (s)
%
% Output:
%       params = parameters for mleBayes, including initial start points 
%       and bounds for each parameter
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

% displacements
deltax = diff(pos);
[N,dim] = size(deltax);

% initial velocity
v0 = mean(deltax)/dt;

% initial diffusivity and static noise estimate (assuming normal diffusion)
C = 0;
for i = 1:dim
    delta = (deltax(:,i)-v0(i)*dt);
    C = C + delta*delta';
end
C = C/dim;
vacf = CalculateVACF(C);
D0 = max((vacf(1) + 2*vacf(2))/2/dt,1e-6);
sigma0 = min(sqrt(abs(vacf(1)-2*D0*dt*(1-2*R))/2),.3);

% initial confinement size estimate
L0 = mean(max(pos) - min(pos));

params.models = models;
params.dt = dt;
params.R = R;
params.D0 = D0;
params.sigma0 = sigma0;
params.v0 = v0;
params.L0 = L0;
params.A0 = 1;

% bounds for parameters
params.Dmin = 1e-8;
params.Dmax = D0*100;
params.Smin = 0;
params.Smax = 1;
params.Vmin = -abs(v0)*10;
params.Vmax = abs(v0)*10;
params.Lmin = L0/10;
params.Lmax = 10;
params.Amin = .1;
params.Amax = 2;

