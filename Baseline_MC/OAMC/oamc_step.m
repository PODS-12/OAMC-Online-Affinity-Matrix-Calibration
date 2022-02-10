function [u_cal] = oamc_step(s, u0, tol)
% function [u_cal] = oamc_step(s, u0, tol)
% 
% @param
% @param

if (nargin < 3)
    tol = 1e-4;
end

if norm(u0) <= 1
    u_cal = u0;
else
    lambda_min = max(sqrt(s'.^2 * u0.^2) - max(s), 0);
    lambda_max = sqrt(s'.^2 * u0.^2) - min(s);
    lambda = lambda_min;
    u_cal = (s ./ (s+lambda)) .* u0;
    u_len = norm(u_cal);
    while (u_len > 1) || (u_len < 1-tol)
        lambda = 0.5*(lambda_min + lambda_max);
        u_cal = (s ./ (s+lambda)) .* u0;
        u_len = norm(u_cal);
        % use bisection search to find optimal lambda
        if u_len > 1
            lambda_min = lambda;
        elseif u_len < 1-tol
            lambda_max = lambda;
        end 
    end
end

end