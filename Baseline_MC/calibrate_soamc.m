function [Soamc] = calibrate_soamc(Smiss, noff, non, model, maxiter)
% function [Soamc] = calibrate_soamc(Smiss, noff, non, model, maxiter)
%
% @param Smiss     Initial similarity matrix
% @param noff      Number of offline samples
% @param non       Number of online samples
% @param koff      size of offline submatrix
% @param kon       size of online submatrix 
% @param model     Type of direct calibration method ('DMC' or 'CMC')
% @param maxiter   Maximum iterations
%
% @return Soamc    Calibrated similarity matrix

if (nargin < 5)
    maxiter = 50;
elseif (nargin < 4)
    model = 'dmc';
end
koff = 1000; kon = 1000;

npart_off = noff / koff;
Soff_par = Smiss(1:noff, noff+1:end);
for i = 1 : npart_off
    scale = [(i-1)*koff+1 : i*koff];
    Soff = Smiss(scale, scale);
    Spar_ini = Smiss(scale, noff+1:end);
    Spar_cal = oamc_parallel(Soff, Spar_ini);
    Soff_par(scale, :) = Spar_cal;
end

Son = Smiss(noff+1:end, noff+1:end);
Son_cal = calibrate_on_para(Son, kon, model, maxiter);

Soamc = [Smiss(1:noff, 1:noff), Soff_par; Soff_par', Son_cal];

end


