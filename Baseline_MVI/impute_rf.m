function [imputedX, time] = impute_rf(Xoff, Xon, model, ntree)
% function [imputedX, time] = impute_rf(Xoff, Xon, model, ntree)
%
% Impute a data matrix. Each column is a sample. Each NaN value in a vector 
% is replaced by the statistical value calculated by Random Forest from
% observed values to missing values. (see reference)
%
% @param Xoff       Offline dataset, each column is a complete sample
% @param Xon        Online dataset, each column is an incomplete sample
% @param model      Default 'on' (online version)
% @param ntree      Default 500 (number of trees)
% 
% @return imputedX  Imputed matrix with all data samples
% @return time      Running time
%
% <Reference>
% Stekhoven, Daniel J., and Peter Bühlmann. "MissForest—non-parametric 
% missing value imputation for mixed-type data." Bioinformatics 28.1 (2012): 112-118.

if (nargin < 4)
    ntree = 500;
elseif (nargin < 3)
    model = 'on';
end
low = -1; high = 1;

n_off = size(Xoff, 2);
n_on = size(Xon, 2);
imputedX = [Xoff, Xon];
if mod(n_on, 100) == 0
    time = zeros(1, n_on/100);
end

if strcmp(model, 'on')
    tic;
    for i = 1 : n_on
        y = Xon(:,i);
        idx = isnan(y);
        y_obs = y(~idx);
        X_mis = Xoff(idx,:);
        X_obs = Xoff(~idx,:);
        RF = TreeBagger(ntree, X_obs, y_obs, 'Method', 'regression');
        y_mis = predict(RF, X_mis);
        y(idx) = y_mis;
        imputedX(:, n_off+i) = y;
        if mod(i, 100) == 0
            time(1, i/100) = toc;
        end
    end
elseif strcmp(model, 'seq')
    tic;
    for i = 1 : n_on
        y = Xon(:,i);
        idx = isnan(y);
        y_obs = y(~idx);
        X_mis = Xoff(idx,:);
        X_obs = Xoff(~idx,:);
        RF = TreeBagger(ntree, X_obs, y_obs, 'Method', 'regression');
        y_mis = predict(RF, X_mis);
        y(idx) = y_mis;
        imputedX(:, n_off+i) = y;
        Xoff = [Xoff, y];
        if mod(i, 100) == 0
            time(1, i/100) = toc;
        end
    end    
elseif strcmp(model, 'off')
    printf('No Offline Version!');
end
imputedX(imputedX<low) = low;
imputedX(imputedX>high) = high;
end