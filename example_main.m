%% 

clear all; clc; warning off;
addpath(genpath('Baseline_MVI'));
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
noff = 1000;               % number of offline samples
non = 100;                % number of online samples
seed = 2022;              % random seed for reproducing results

fprintf('\nSIGKDD 2022 submission #730 "Online Affinity Matrix Calibration with Incomplete Observations"');
fprintf('\nDemo: missing data processing in Section 5.3\n');

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
        
        %% Online Process
        % ====================== ZERO Imputation ==========================
        fprintf('ZERO, ');
        tic; Xzero = impute_zero(Xmiss); time(i,1) = toc;
        Szero = similarity(Xzero);
        rmse(i,1) = norm(Szero-Strue, 'fro')^2 / Fnorm;
        
        % ====================== MEAN Imputation ==========================
        fprintf('MEAN, ');
        [Xmean, Tmean] = impute_mean(Xoff, Xon, 'on'); time(i,2) = Tmean(end);
        Smean = similarity(Xmean);
        rmse(i,2) = norm(Smean-Strue, 'fro')^2 / Fnorm;

        % ====================== kNN Imputation ===========================
        fprintf('kNN, '); 
        [Xknn, Tknn] = impute_knn(Xoff, Xon, 'on'); time(i,3) = Tknn(end);
        Sknn = similarity(Xknn);
        rmse(i,3) = norm(Sknn-Strue, 'fro')^2 / Fnorm;
        
        % ====================== LR Imputation ============================
        fprintf('LR, '); 
        [Xlr, Tlr] = impute_lr(Xoff, Xon, 'on'); time(i,4) = Tlr(end);
        Slr = similarity(Xlr);
        rmse(i,4) = norm(Slr-Strue, 'fro')^2 / Fnorm;
        
        % ====================== RF Imputation ============================
        fprintf('RF, '); 
        [Xrf, Trf] = impute_rf(Xoff, Xon, 'on', 50); time(i,5) = Trf(end);
        Srf = similarity(Xrf);
        rmse(i,5) = norm(Srf-Strue, 'fro')^2 / Fnorm;
        
        % ====================== GROUSE Imputation ========================
        fprintf('GROUSE, '); 
        [Xgr, Tgr] = impute_grouse(Xoff, Xon, 'off'); time(i,6) = Tgr(end);
        Sgr = similarity(Xgr);
        rmse(i,6) = norm(Sgr-Strue, 'fro')^2 / Fnorm;
        
        % ====================== KFMC Imputation ==========================
        fprintf('KFMC, '); 
        tic; Xkfmc = impute_kfmc(Xmiss, Xtrue, 'on', 'rbf'); time(i,7) = toc;
        Skfmc = similarity(Xkfmc);
        rmse(i,7) = norm(Skfmc-Strue, 'fro')^2 / Fnorm;
        
        % ====================== DMC Calibration ==========================
        fprintf('DMC, '); 
        tic; Sdmc = calibrate_dmc(Smiss); time(i,8) = toc;
        rmse(i,8) = norm(Sdmc-Strue, 'fro')^2 / Fnorm;

        % ====================== CMC Calibration ==========================
        fprintf('CMC, '); 
        tic; Scmc = calibrate_cmc(Smiss); time(i,9) = toc;
        rmse(i,9) = norm(Scmc-Strue, 'fro')^2 / Fnorm;
        
        % ====================== OAMC Calibration =========================
        fprintf('OAMC, '); 
        [Soamc, Toamc] = calibrate_oamc(Smiss, noff, non); time(i,10) = Toamc(end);
        rmse(i,10) = norm(Soamc-Strue, 'fro')^2 / Fnorm;
        
        % ====================== OAMC-DMC Calibration =====================
        fprintf('OAMC_DMC, '); 
        tic; Soamc = calibrate_oamc_block(Smiss, noff, non, 'dmc'); time(i,11) = toc;
        rmse(i,11) = norm(Soamc-Strue, 'fro')^2 / Fnorm;
        
        % ====================== OAMC-CMC Calibration =====================
        fprintf('OAMC_CMC, '); 
        tic; Soamc = calibrate_oamc_block(Smiss, noff, non, 'cmc'); time(i,12) = toc;
        rmse(i,12) = norm(Soamc-Strue, 'fro')^2 / Fnorm;
        
        fprintf('Finish.');
    end

    %%
    fprintf(['\n',dataset]); fprintf(': noff=%1.0f, non=%1.0f, p=%1.2f, niter=%1.0f\n', noff, non, p, niter);

    stat = [mean(rmse); mean(time)];
    Stat = roundn(stat, -3);
    Table = table(Stat(:,1),Stat(:,2),Stat(:,3),Stat(:,4),Stat(:,5),Stat(:,6),Stat(:,7),Stat(:,8),Stat(:,9),Stat(:,10),Stat(:,11),Stat(:,12),...
         'VariableNames',{'ZERO','MEAN','kNN','LR','RF','GROUSE','KFMC','DMC','CMC','OAMC','OAMC_DMC','OAMC_CMC'},...
         'RowNames',{'RMSE';'Time'});
    disp(Table)
end



