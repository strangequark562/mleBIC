function  logL = LogMultiVariateNormal(N,deltax,C,mu)
%-------------------------------------------------------------------------- 
% Summary: This function calculates log of a multivariate Normal distribution
% for a 2D particle trajectory.
%
% Input:
%         N = number of displacements
%         deltax = vector of x displacements
%         deltay = vector of y displacements
%         mu = vector of drift velocity*dt
%         C = matrix of covariance matrix 
%
% Output:
%       logL = log-likelihood
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

if nargin == 3
    mu = zeros(1,size(deltax,2));
end

try
    [invC,logdetC] = MatrixInverse(C);
    
    logL = 0;
    for i = 1:size(deltax,2)
        deltax2 = deltax(:,i)-mu(i);
        logL = logL - N/2*log(2*pi) -.5*logdetC - .5*deltax2'*invC*deltax2;
    end
catch
    logL = -Inf;
end

