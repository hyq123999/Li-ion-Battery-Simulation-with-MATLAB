function [ params ] = getParameters()

% Faraday Constant  [C/mol]
params.F = 96487;
% gas constant  [J/mol/K]
params.R = 8.314;

% temperature [T] % assume constant...
params.T = 298.15;

params.alpha = 0.5;

% Legnth of anode/separator/cathode
params.len_c = 80e-6; % [m]
params.len_s = 25e-6; % [m]
params.len_a = 90e-6; % [m]

% Current density at 1C rate
params.i_1C = 29.23; % [A/m^2]

% Solid conductivities [S/m]
params.sigma_c = 100;
params.sigma_a = 100;

% Electrolyte diffusion coefficient
params.De = 7.5e-10; % [m^2/s]

% Electrolyte porosity
params.eps_c = 0.385;
params.eps_s = 0.724;
params.eps_a = 0.485;

% Bruggeman coefficient
params.brug_c = 4.0;
params.brug_s = 4.0;
params.brug_a = 4.0;

% (Li ion) transference number
params.tran_num_plus = 0.364;

% Solid diffusion coefficients 
params.Ds_c = 1.0e-14; % [m^2/s]
params.Ds_a = 3.9e-14; % [m^2/s]

% Particle surface area
params.as_c = 8.85e5; % [m^2/m^3]
params.as_a = 7.236e5; % [m^2/m^3]

% Reaction rate constants [ m ^ 2.5 / (mol^0.5 s ) ]
params.k_c = 2.334e-11; 
params.k_a = 5.031e-11; 

% Maximum concentration of Li-ions in the solid phase [ mol/m^3 ]
params.cs_c_max   = 5.1554e4;
params.cs_a_max   = 3.0555e4;

% Solid particle radius [m] 
params.Rs_c = 2.e-6;
params.Rs_a = 2.e-6;


end




