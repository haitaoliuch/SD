 function mu = SimMeasure(beta, betaHat,problem)
% beta is the true parameter, and betaHat is the estimate obtained by MLE
% Similarity measure is defined as the the size of the intersection divided ...
% by the size of the union of the sample sets with equal weights
% clc,clear,close all
%% ---- Method 2 starts: vectorization to speed up checking feasibility 
%---Define the logistic regression function

% f=@(x,beta) 1./(1+exp(-x' * beta));

X = [ones(size(problem.DS,1),1) problem.DS];
Ytrue = problem.fun(X',beta);
Ypred = problem.fun(X',betaHat);

TP = sum(Ytrue >= problem.eta & Ypred >= problem.eta);
exceed = sum(Ytrue < problem.eta  & Ypred <  problem.eta); % T P + F P + F N

% Compute the similarity measure
if (size(problem.DS,1)-exceed )~=0
    mu= TP/(size(problem.DS,1)-exceed);
else
    mu=0;
end
end
