function [logdetC]= LogDeterminant(C)
%-------------------------------------------------------------------------- 
% Summary: This function calculates log determinant using eigenvalue 
% decomposition to get around numerical underflow issues.
% 
% Input:
%       C = covariance matrix
%
% Output:
%       logdetC = log-determinant
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

threshold = 1e-10;
[vec val] = eig(C);
val = diag(val);
index = find(val <=threshold);
val(index) = threshold;
logdetC = sum(log(val));
