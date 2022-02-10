function [S] = calibrate_dmc(S0, maxiter)
% function [S] = calibrate_dmc(S0, maxiter)
%
% Calibrate a similarity matrix by a direct calibration method (see reference)
%
% @param  S0   Initial similarity matrix
% @return S    Calibrated matrix
%
% <Reference>
% Li, Wenye. "Estimating Jaccard index with missing observations: a matrix 
% calibration approach." Advances in Neural Information Processing Systems 28 (2015).

if (nargin < 2)
    maxiter = 20;
end

low = -1; high = 1;
S = nearpsd(S0, maxiter, low, high);

end

%%
function [X, iter] = nearpsd(A, maxits, low, high)
% function [X, iter] = nearpsd(A, maxits, low, high)
%
% Computes the nearest positive semi-definite matrix 
% for a given square matrix.
%
% @param A        a square matrix to be calibrated
% @param maxits   max num of iters allowed, default 100
% @param low      default 0 
% @param high     default 1
%
% @return X       nearest psd matrix to A
% @return iter    number of iterations taken
%

if  ~isequal(A,A')
    A = (A + A') / 2;
end
if nargin < 4
    high = 1;
end
if nargin < 3
    low = -1; 
end
if nargin < 2
    maxits = 20;
end

% threshold for convergence & eigs
tolconv = 1.0e-6;
toleigs = 1.0e-5;

n = size(A,1);

U = zeros(n);
Y = A;

[V, D] = eig(Y);
d = diag(D);

iter = 0;
while 1
    T = Y - U;

    % project onto psd matrices
    [Q, D] = eig(T);
    d = diag(D);
    p = d > toleigs*d(n);
    X = Q(:,p) * D(p,p) * Q(:,p)';

    % update correction
    U = X - T;

    % maximum iteration & convergence test
    iter = iter + 1;
    if iter == maxits
        %fprintf('Max iterations reached. ');
        break; 
    end
    if norm(Y-X,'inf')/norm(Y,'inf') <= tolconv 
    	break;
    end
    
    % problem-dependent here
    Y = X;
    Y(1:n+1:n*n) = 1;
    Y(Y<low) = low;
    Y(Y>high) = high;
end

Y(1:n+1:n*n) = 1;
Y(Y<low) = low;
Y(Y>high) = high;
end
