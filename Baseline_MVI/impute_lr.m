function [imputedX, time] = impute_lr(Xoff, Xon, model)
% function [imputedX, time] = impute_lr(Xoff, Xon, model)
%
% Impute a data matrix. Each column is a sample. Each NaN value in a vector 
% is replaced by the statistical value calculated by Linear Regression from
% observed values to missing values. (see reference)
%
% @param Xoff       Offline dataset, each column is a complete sample
% @param Xon        Online dataset, each column is an incomplete sample
% @param model      Default 'on' (online version)
% 
% @return imputedX  Imputed matrix with all data samples
% @return time      Running time
%
% <Reference>
% Seber, George AF, and Alan J. Lee. Linear regression analysis. John Wiley & Sons, 2012.

if (nargin < 3)
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
        beta = regress(y_obs, [ones(sum(~idx),1), X_obs]);
        y_mis = sum(beta' .* [ones(sum(idx),1), X_mis], 2);
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
        beta = regress(y_obs, [ones(sum(~idx),1), X_obs]);
        y_mis = sum(beta' .* [ones(sum(idx),1), X_mis], 2);
        y(idx) = y_mis;
        imputedX(:, n_off+i) = y;
        Xoff = [Xoff, y];
        if mod(i, 100) == 0
            time(1, i/100) = toc;
        end
    end    
elseif strcmp(model, 'off')
    fprintf('No Offline Version!');
end
imputedX(imputedX<low) = low;
imputedX(imputedX>high) = high;
end