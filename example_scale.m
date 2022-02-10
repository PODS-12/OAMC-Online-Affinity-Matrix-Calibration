%% 

clear all; clc;
addpath(genpath('Baseline_MC'));

load('demo_data.mat');
dataset = 'Demo_data'; 

%% Data Normalization
X = double(X);
n = size(X,2);
% normalize to [-1, 1]
X = X - repmat(min(X,[],2), 1, n);
X = X ./ (repmat(max(X,[],2), 1, n) + eps);
X = (X - 0.5) * 2;

%% Parameter Setting
p_list = [0.2, 0.5, 0.8]; % values of missing ratio p
niter = 5;                % number of iterations
noff = 5000;              % number of offline samples
non = 5000;               % number of online samples
seed = 2022;              % random seed for reproducing results

fprintf('\nSIGKDD 2022 submission #730 "Online Affinity Matrix Calibration with Incomplete Observations"');
fprintf('\nDemo: scalability analysis in Section 5.4\n');

for k = 1 : length(p_list)
    p = p_list(k);
    fprintf(['\n',dataset]); fprintf(': noff=%1.0f, non=%1.0f, p=%1.2f, niter=%1.0f', noff, non, p, niter);
    rng(seed);
    
    %% Online Scenario
    for i = 1 : niter
        fprintf('\nIters = %1.0f: ', i);
        
        %% Data Construction
        rp = randperm(n);
        Xoff = X(:, rp(1:noff));
        Xon = X(:, rp(noff+1 : noff+non));
        Xtrue = [Xoff, Xon];
        Strue = similarity(Xtrue);

        Xon(rand(size(Xon)) < p) = NaN;  
        Xmiss = [Xoff, Xon];
        Smiss = similarity(Xmiss, 'miss');
        Fnorm = norm(Smiss-Strue, 'fro')^2;
        
        %% SOAMC Calibration
        % ====================== SOAMC-DMC Calibration ====================
        fprintf('SOAMC_DMC, '); 
        tic; Soamc = calibrate_soamc(Smiss, noff, non, 'dmc', 10); time(i,1) = toc;
        rmse(i,1) = norm(Soamc-Strue, 'fro')^2 / Fnorm;
        
        % ====================== SOAMC-CMC Calibration ====================
        fprintf('SOAMC_CMC, '); 
        tic; Soamc = calibrate_soamc(Smiss, noff, non, 'cmc', 10); time(i,2) = toc;
        rmse(i,2) = norm(Soamc-Strue, 'fro')^2 / Fnorm;
        
        fprintf('Finish.');
    end

    %%
    fprintf(['\n',dataset]); fprintf(': noff=%1.0f, non=%1.0f, p=%1.2f, niter=%1.0f\n', noff, non, p, niter);

    stat = [mean(rmse); mean(time)];
    Stat = roundn(stat, -3);
    Table = table(Stat(:,1),Stat(:,2),...
         'VariableNames',{'SOAMC_DMC','SOAMC_CMC'},...
         'RowNames',{'RMSE';'Time'});
    disp(Table)
end



