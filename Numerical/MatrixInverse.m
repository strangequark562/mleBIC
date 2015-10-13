function [invC,logdetC] = MatrixInverse(C)
%-------------------------------------------------------------------------- 
% Summary: This function calculates matrix inverse and log determinant.
% This function uses eigenvalue decomposition to get around numerical
% underflow issues.
% 
% Input:
%       C = covariance matrix
%
% Output:
%       invC = inverse covariance matrix
%       logdetC = log-determinant
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

threshold = 1e-10;
[vec,val] = eig(C);
val = diag(val);
index = find(val <=threshold);
val(index) = threshold;
logdetC = sum(log(val));
invC = vec*diag(1./val)*vec';
    
    
