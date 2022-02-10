function [S] = calibrate_on_para(Son, kon, model, maxiter)
% function [S] = calibrate_on_para(Son, kon, model, maxiter)
%
% @param Son       Initial online matrix
% @param kon       size of online submatrix 
% @param model     Type of direct calibration method ('DMC' or 'CMC')
% @param maxiter   Maximum iterations
%
% @return S        Calibrated online matrix

if (nargin < 4)
    maxiter = 50;
elseif (nargin < 3)
    model = 'dmc';
elseif (nargin < 2)
    kon = 1000;
end

non = size(Son, 1);
npart_on = non / kon;
S = zeros(non, non);
for i = 1 : npart_on
    scale = [(i-1)*kon+1 : i*kon];
    Son_ini = Son(scale, scale);
    Spar_ini = Son(scale, i*kon+1:non);
    if strcmp(model, 'dmc')
        Son_cal = calibrate_dmc(Son_ini, maxiter);
    elseif strcmp(model, 'cmc')
        Son_cal = calibrate_cmc(Son_ini, maxiter);
    end
    Spar_cal = oamc_parallel(Son_cal, Spar_ini);
    S(scale, (i-1)*kon+1:end) = [Son_ini, Spar_cal];
    S(i*kon+1:end, scale) = Spar_cal';
end

end


