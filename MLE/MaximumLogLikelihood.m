function  [mu,I,maxLogL] = MaximumLogLikelihood(pos,params)
%-------------------------------------------------------------------------- 
% Summary: This function finds the maximum log-likelihood estimates for the 
% diffusion models given by params.models
% 
% Input:
%       x = x positions of a single particle trajectory in 2D
%       y = y positions of a single particle trajectory in 2D
%       params = parameters for MLE fitting (see InitialParameters.m)
%
% Output:
%       mu = cell of parameter vectors of each model 
%       I = cell of empirical Fisher information matrix for each model 
%       maxLogL = vector of log-likelihood values with elements 
%                 correpsonding to  each diffusion model
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

% optimization parameters 'algorithm','sqp',
fitOptions = optimset('Display','off','TolX',1e-6,'MaxFunEvals',1000,'MaxIter',1000);

% particle track displacements
deltax = diff(pos);
[N,dim] = size(deltax);

% initial parameters
models = params.models;
dt = params.dt;
R = params.R;
D0 = params.D0;
sigma0 = params.sigma0;
v0 = params.v0;
L0 = params.L0;
A0 = params.A0;

% bounds for parameters
Dmin = params.Dmin;
Dmax = params.Dmax;
Smin = params.Smin;
Smax = params.Smax;
Vmin = params.Vmin;
Vmax = params.Vmax;
Lmin = params.Lmin;
Lmax = params.Lmax;
Amin = params.Amin;
Amax = params.Amax;

% loop through each diffusion model and find MLE estimates
numModels = length(models);
mu = cell(numModels,1); 
I = cell(numModels,1);
maxLogL = zeros(numModels,1);
for j = 1:numModels
    switch models{j}
        case 'Immobile'
            lb = [Smin];
            ub = [Smax];
            initialP0 = [sigma0];
            [p,logL,e,o,l,g,hessian] = fmincon(@(b)-Likelihood_Immobile(N,deltax,b),initialP0,[],[],[],[],lb,ub,[],fitOptions);
        
        case 'Normal'
            lb = [Dmin Smin];
            ub = [Dmax Smax];
            initialP0 = [D0 sigma0];
            [p,logL,e,o,l,g,hessian] = fmincon(@(b)-Likelihood_NormalDiffusion(N,deltax,b(1),b(2),R,dt),initialP0,[],[],[],[],lb,ub,[],fitOptions);

        case 'Driven'
            lb = Dmin;
            ub = Dmax;
            for i = 1:dim
                lb = [lb Vmin];
                ub = [ub Vmax];
            end
            lb = [lb Smin];
            ub = [ub Smax];
            initialP0 = [D0/2 v0 sigma0];
            [p,logL,e,o,l,g,hessian] = fmincon(@(b)-Likelihood_DrivenDiffusion(N,deltax,b,R,dt),initialP0,[],[],[],[],lb,ub,[],fitOptions);

        case 'Confined'
            lb = [Dmin Lmin Smin];
            ub = [Dmax Lmax Smax];
            if L0 < .5
                initialP0 = [D0*2 L0 sigma0];
            else
                initialP0 = [D0 L0 sigma0];
            end
            [p,logL,e,o,l,g,hessian] = fmincon(@(b)-Likelihood_confinedDiffusion(N,deltax,b(1),b(2),b(3),R,dt),initialP0,[],[],[],[],lb,ub,[],fitOptions);

        case 'fBM'
            lb = [Dmin Amin Smin];
            ub = [Dmax Amax Smax];
            if sigma0 > .1
                initialP0 = [D0/2 A0*.6 .04];
            elseif sigma0 < .01
                initialP0 = [D0*2 A0*1.3 sigma0*2];
            end
            [p,logL,e,o,l,g,hessian] = fmincon(@(b)-Likelihood_fBM(N,deltax,b(1),b(2),b(3),R,dt),initialP0,[],[],[],[],lb,ub,[],fitOptions);
        
        otherwise
            disp('no case found');
    end
    
    % store estimates if doesn't return junk
    mu{j} = p';
    if isempty(logL)
        maxLogL(j) = -Inf;
        I{j} = 1;
    else
        maxLogL(j) = -logL;
        I{j} = MatrixInverse(hessian);
    end
end


