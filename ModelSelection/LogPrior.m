function logPrior = LogPrior(p,mu,C)
%-------------------------------------------------------------------------- 
% Summary: This function calculates log priors of the parameters
% 
% Input:
%       p = parameters
%       mu = hyperparameter of the multivariate Gaussian
%       C = hyperparameter of the multivariate Gaussian
%
% Output:
%       logPrior = log prior value
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

if isempty(C)
    logPrior = -Inf;
else
    [invC,logdetC] = MatrixInverse(C);    
    deltax2 = (p-mu);
    logPrior = -length(p)/2*log(2*pi) -.5*logdetC - .5*deltax2'*invC*deltax2;
end
