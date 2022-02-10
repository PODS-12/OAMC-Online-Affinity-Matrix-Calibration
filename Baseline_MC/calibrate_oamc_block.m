function [Soamc] = calibrate_oamc_block(Smiss, noff, non, model, maxiter)
% function [Soamc] = calibrate_oamc_block(Smiss, noff, non, model, maxiter)
%
% @param Smiss     Initial similarity matrix
% @param noff      Number of offline samples
% @param non       Number of online samples
% @param model     Type of direct calibration method ('DMC' or 'CMC')
% @param maxiter   Maximum iterations
%
% @return Soamc    Calibrated similarity matrix

if (nargin < 5)
    maxiter = 50;
elseif (nargin < 4)
    model = 'cmc';
end

Son = Smiss(noff+1:end, noff+1:end);
if strcmp(model, 'dmc')
    Son_cal = calibrate_dmc(Son, maxiter);
elseif strcmp(model, 'cmc')
    Son_cal = calibrate_cmc(Son, maxiter);
end

Soff = Smiss(1:noff, 1:noff);
Spar_ini = Smiss(1:noff, noff+1:end);
Spar_cal = oamc_parallel(Soff, Spar_ini);

Soamc = [Smiss(1:noff, 1:noff), Spar_cal; Spar_cal', Son_cal];

end


