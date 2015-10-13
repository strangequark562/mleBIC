function  [mu,I,maxLogL] = MaximumLogLikelihoodNoNoise(x,y,params)

% optimization parameters 'algorithm','sqp',
fitOptions = optimset('Display','off','TolX',1e-6,'MaxFunEvals',1000,'MaxIter',1000);

% particle track displacements
deltax = diff(x);
deltay = diff(y);
N = length(deltax);

% initial parameters
models = params.models;
dt = params.dt;
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
            [p,logL,e,o,l,g,hessian] = fmincon(@(b)-Likelihood_Immobile(N,deltax,deltay,b),initialP0,[],[],[],[],lb,ub,[],fitOptions);
        case 'Normal'
            lb = [Dmin];
            ub = [Dmax];
            initialP0 = [D0];
            [p,logL,e,o,l,g,hessian] = fmincon(@(b)-Likelihood_NormalDiffusion(N,deltax,deltay,b(1),0,0,dt),initialP0,[],[],[],[],lb,ub,[],fitOptions);

        case 'Driven'
            lb = [Dmin Vmin];
            ub = [Dmax Vmax];
            initialP0 = [D0/2 v0];
            [p,logL,e,o,l,g,hessian] = fmincon(@(b)-Likelihood_DrivenDiffusion(N,deltax2,deltay2,b(1),b(2),0,0,dt),initialP0,[],[],[],[],lb,ub,[],fitOptions);
        
        case 'Confined'
            lb = [Dmin Lmin];
            ub = [Dmax Lmax];
            if L0 < .5
                initialP0 = [D0*2 L0];
            else
                initialP0 = [D0 L0];
            end
            [p,logL,e,o,l,g,hessian] = fmincon(@(b)-Likelihood_confinedDiffusion(N,deltax,deltay,b(1),b(2),0,0,dt),initialP0,[],[],[],[],lb,ub,[],fitOptions);
        
        case 'fBM'
            lb = [Dmin Amin];
            ub = [Dmax Amax];
            if sigma0 > .1
                initialP0 = [D0/2 A0*.7];
            elseif sigma0 < .01
                initialP0 = [D0*2 A0*1.3];
            end
            [p,logL,e,o,l,g,hessian] = fmincon(@(b)-Likelihood_fBM(N,deltax,deltay,b(1),b(2),0,0,dt),initialP0,[],[],[],[],lb,ub,[],fitOptions);
        
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


