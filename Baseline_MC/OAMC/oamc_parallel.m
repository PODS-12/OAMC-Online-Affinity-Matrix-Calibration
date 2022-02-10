function [Son_cal] = oamc_parallel(Soff, Son, s, V)
% function [Son_cal] = oamc_parallel(Soff, Son)
%
% @param  S0    pairwise similarity matrix
% @param  v0    intial similarity vector
% @return vonl  calibrated similarity vector

if (nargin < 4)
    [V, S, ~] = svd(Soff);
    s = diag(S);
end

n = size(Son, 2);
tol = 1e-4;
C = V * diag(sqrt(s));
Cinv = diag(1./s) * C';
U = Cinv * Son;
U_cal = U;

for i = 1 : n
    u0 = U(:, i);
    U_cal(:, i) = oamc_step(s, u0, tol);
end
Son_cal = C * U_cal;

end

