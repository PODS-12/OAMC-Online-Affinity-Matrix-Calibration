function [imputedX, time] = impute_knn(Xoff, Xon, model, k)
% function [imputedX, time] = impute_knn(Xoff, Xon, k, model)
%
% Impute a data matrix. Each column is a sample. Each NaN value is replaced
% by the mean of the sample's k-nearest neighbors with known values. If all 
% k-nearest samples' corresponding features are NaN, then replaced by zero.
%
% @param Xoff       Offline dataset, each column is a complete sample
% @param Xon        Online dataset, each column is an incomplete sample
% @param model      Default 'on' (online version)
% @param k          Default 10 (k-nearest neighbors)
% 
% @return imputedX  Imputed matrix with all data samples
% @return time      Running time

if (nargin < 4)
    k = 10;
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
        X = [Xoff, Xon(:, i)];
        Ximp = knnimpute(X, k);
        imputedX(:, n_off+i) = Ximp(:, end);
        if mod(i, 100) == 0
            time(1, i/100) = toc;
        end
    end
elseif strcmp(model, 'seq')
    tic;
    for i = 1 : n_on
        X = [Xoff, Xon(:, i)];
        Ximp = knnimpute(X, k);
        imputedX(:, n_off+i) = Ximp(:, end);
        Xoff = [Xoff, Ximp(:, end)];
        if mod(i, 100) == 0
            time(1, i/100) = toc;
        end
    end    
elseif strcmp(model, 'off')
    tic;
    imputedX = SeqKNN(imputedX', k)';
    time = toc;
end
imputedX(imputedX<low) = low;
imputedX(imputedX>high) = high;
end


%%
function [imputed_data] = SeqKNN(data, K)
% SeqKNN: Sequential KNN imputation method
% This function estimates missing values sequentially from the gene that has
% least missing rate in microarray data, using weighted mean of k nearest neighbors.
%
% <Usage>
% imputed_data = SeqKNN(data, k);
%
% <Arguments>
% data: matrix or dataframe, 1 row corresponds to 1 gene, 1 column to 1
% sample,colnames and rownames can be used
% k: number of nearest neighbors
%
% <Details>
% SeqKNN� separates the dataset into incomplete and complete set that has or has not missing values
% respectively. The genes in incomplete set are imputed by the order of missing rate. Missing value
% is filled by the weighted mean value of corresponding column of the nearest neighbor genes in
% complete set. Once all missing values in a gene are imputed, the imputed gene is moved into the
% complete set and used for the imputation of the rest of genes in incomplete set. In this process,
% all missing values in one gene can be imputed simultaneously from the selected neighbor genes
% in complete set. This reduces execution time from previously developed KNN method that selects
% nearest neighbors for each imputation.
%
% <Reference>
% Ki-Yeol Kim, Byoung-Jin Kim, Gwan-Su Yi (2004.Oct.26) "Reuse of imputed data in microarray
% analysis increases imputation efficiency", BMC Bioinformatics 5:160.
% -------------------------------------------------------------------------
% Part 1
% display(sprintf('Working ......'));
imputed_data=zeros(size(data));
complete=[];
incomplete=[];
missing=[];
com_ind=[];
incom_ind=[];
[rows,cols]=size(data);
for i=1:rows
    if ~isnan(sum(data(i,:)))
        complete=[complete; data(i,:)];
        com_ind=[com_ind i];
    else
        incomplete=[incomplete; data(i,:)];
        incom_ind=[incom_ind i];
        missing=[missing sum(isnan(data(i,:)))];
    end
end
imputed_data(com_ind,:)=complete;

[missing,missing_order]=sort(missing);
% incomplete=incomplete(missing_order,:);

% Part2
% K=10;
[irows,icols]=size(incomplete);
for j=1:irows
    dist=[];
    cgen=size(complete,1);
    for i=1:cgen
        dist(i)=nansum((incomplete(missing_order(j),:)-complete(i,:)).^2);
    end
    [dist,pos]=sort(dist);
    pos=pos(1:K);
    dist=dist(1:K);
    weights=(1./dist)./sum(1./dist);
    for g=1:icols
        if (isnan(incomplete(missing_order(j),g)))
%             incomplete(missing_order(j),g)=mean(complete(pos,g));
            incomplete(missing_order(j),g)=weights*(complete(pos,g));
        end
    end
    complete=[complete;incomplete(missing_order(j),:)];
end
imputed_data(incom_ind,:)=incomplete;
% display(sprintf('------ DONE ------'));
end