function [imputedX] = impute_kfmc(Xmiss, Xtrue, model, method)
% function [imputedX] = impute_kfmc(Xmiss, Xtrue, model, method)
%
% Impute a data matrix. Each column is a sample. Each NaN value is replaced
% by the statistical value calculated by the KFMC method (see reference).
%
% @param Xmiss     All original data including offline data and online data
% @param Xtrue     The ground truth of all data samples
% @param model     Default 'on' (online version)
% @param method    Default 'rbf' (RBF kernel function)
% 
% @return imputedX Imputed matrix with all data samples
%
% <Reference>
% Fan, Jicong, and Madeleine Udell. "Online high rank matrix completion." 
% Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition. 2019.

if (nargin < 4)
    method = 'rbf';
elseif (nargin < 3)
    model = 'on';
end
low = -1; high = 1;

Xzero = Xmiss;
Xzero(isnan(Xmiss)==1) = 0;
M = double(~isnan(Xmiss));

if strcmp(model, 'on')
    if strcmp(method, 'poly')
        ker.type = 'poly'; ker.par = [1, 2]; d = 100;
        options.online_maxiter = 20; options.eta = 0.5; options.X_true = Xtrue; options.npass = 2; 
        [imputedX, ~] = KFMC_online(Xzero, M, d, 0.01, 0.001, ker, options);
    elseif strcmp(method, 'rbf')
        ker.type = 'rbf'; ker.par = 0; ker.c = 3; d = 100;
        options.online_maxiter = 20; options.eta = 0.5; options.X_true = Xtrue; options.npass = 2;
        [imputedX, ~] = KFMC_online(Xzero, M, d, 0, 0.001, ker, options);
    end
elseif strcmp(model, 'off')
    if strcmp(method, 'poly')
        ker.type = 'poly'; ker.par = [1 2];
        alpha = 0; beta = 0.001; d = 100;
        [imputedX, ~] = KFMC(Xzero, M, d, alpha, beta, ker);
    elseif strcmp(method, 'rbf')
        ker.type = 'rbf'; ker.par = []; ker.par_c = 1;
        alpha = 0; beta = 0.001; d = 100;
        [imputedX, ~] = KFMC(Xzero, M, d, alpha, beta, ker);
    end
end  
imputedX(imputedX<low) = low;
imputedX(imputedX>high) = high;
end