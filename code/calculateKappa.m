function [ kappa ] = calculateKappa(ce, params)
% % ce is a vector
% % kappa is a vector
    kappa = params.eps .^ params.brug .*(4.1253*1e-2 + 5.007*1e-4*ce - 4.7212*1e-7*ce.^2 +1.5094*1e-10*ce.^3 -1.6018*1e-14*ce.^4);
end