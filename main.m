%-------------------------------------------------------------------------- 
% Summary: This script prompts the user to select a .mat file which
% contains a cell, X, which contains each particle trajectory:
% X{1} = [x1,y1]; X{2} = [x2,y2]; ... X{n} =[xn,yn], where xi and yi are
% vectors positions for particle track i. The output is a contained in
% results.  A summary of the results is displayed on the command prompt.
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 
clear all;
clc;
close all;
warning off;

addpath('MLE');
addpath('ModelSelection');
addpath('Numerical');

%%

% load files
[filename,dirpath] = uigetfile('*.mat','select file');
load(fullfile(dirpath,filename));

% model parameters
models = {'Immobile','Normal','Driven','Confined','fBM'};
R = 1/6;

% perform mleBayes on each trajectory
k = 1;
results = struct([]);
for n = 1:length(X)
    % mleBayes model selection
%     strength = 1e3;
%     results(n) = mleBayes(X{n},models,R,dt,strength);
    
    % BIC model selection
    results(n) = mleBIC(X{n},models,R,dt);
end



%% 



