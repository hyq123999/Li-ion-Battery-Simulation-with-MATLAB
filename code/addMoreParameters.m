function [ params ] = addMoreParameters( params, rate, h, hr, k )
%% 
params.IoverA = params.i_1C * rate; % scalar

params.h = h;
params.hr = hr;
params.k = k; 

%%
params.m_c = round(params.len_c/params.h) - 1; 
params.m_s = round(params.len_s/params.h) - 1; 
params.m_a = round(params.len_a/params.h) - 1;
params.m = (params.m_c+params.m_s+params.m_a) + 2; % number of "interior" spatial points 

params.mr_c = (params.Rs_c*1e8)/(params.hr*1e8) - 1;
params.mr_a = (params.Rs_a*1e8)/(params.hr*1e8) - 1;

params.eps = [ones(params.m_c+2,1)*params.eps_c;...
              ones(params.m_s+0,1)*params.eps_s;...
              ones(params.m_a+2,1)*params.eps_a];
params.brug = [ones(params.m_c+2,1)*params.brug_c;...
               ones(params.m_s+0,1)*params.brug_s;...
               ones(params.m_a+2,1)*params.brug_a];

%% electrolyte           
De = params.De;

m_c = params.m_c; m_s = params.m_s; m_a = params.m_a;
m = params.m; % number of "interior" spatial points 
 
eps = params.eps;
brug = params.brug;

De_eff = De * eps .^ brug;

% sharp jump??
% % a = diag([-3; 1/k*eps(2:(end-1)) + 2/h^2*De_eff(2:(end-1)); -3])...
% %     +diag([4; -1/h^2*De_eff(3:end)],1)...
% %     +diag([-1/h^2*De_eff(1:(end-2)); 4],-1);
a = diag([-3; 1/k*eps(2:(end-1)) + 2/h^2*De_eff(2:(end-1)); -3])...
    +diag([4; -1/h^2*De_eff(2:(end-1))],1)...
    +diag([-1/h^2*De_eff(2:(end-1)); 4],-1);
%
a(1,3) = -1; 
a(m+2,m) = -1;
%
ib1 = m_c+2;
a(ib1,ib1-2) = -1*De_eff(ib1-1); a(ib1,ib1-1) =  4*De_eff(ib1-1); 
a(ib1,ib1) = -3*De_eff(ib1-1) -3*De_eff(ib1+1);
a(ib1,ib1+1) =  4*De_eff(ib1+1); a(ib1,ib1+2) = -1*De_eff(ib1+1); 
%
ib2 = m_c+m_s+3;
a(ib2,ib2-2) = -1*De_eff(ib2-1); a(ib2,ib2-1) =  4*De_eff(ib2-1); 
a(ib2,ib2) = -3*De_eff(ib2-1) -3*De_eff(ib2+1);
a(ib2,ib2+1) =  4*De_eff(ib2+1); a(ib2,ib2+2) = -1*De_eff(ib2+1); 

b = zeros(m+2,1); 
S = zeros(m+2,1);

params.ib1 = ib1;
params.ib2 = ib2;
params.speciesAtElectrolyte_a = a;
params.speciesAtElectrolyte_b = b;
params.speciesAtElectrolyte_S = S;
           
%% cathode (of electordes)
Ds = params.Ds_c; Rs = params.Rs_c;
mr = params.mr_c;

sigma = Ds * k / (hr*hr); 
ri = (hr:hr:(Rs-hr))'; 
ri_minus_half = -1/2*hr + ri; ri_plus_half = 1/2*hr + ri;
Ri_minus_square = (ri_minus_half./ri).^2; Ri_plus_square = (ri_plus_half./ri).^2;

a=diag([-3; 1 + sigma*(Ri_minus_square+Ri_plus_square); -3])...
    +diag([4; -sigma*Ri_plus_square], 1)...
    +diag([-sigma*Ri_minus_square; 4], -1);
a(1,3) = -1; a(mr+2,mr) = -1;

b = zeros(mr+2,1); S = zeros(mr+2,1);

params.speciesAtSolid_Cathode_a = a;
params.speciesAtSolid_Cathode_b = b;
params.speciesAtSolid_Cathode_S = S;

%% anode (of electrodes)
Ds = params.Ds_a; Rs = params.Rs_a;
mr = params.mr_a;

sigma = Ds * k / (hr*hr); 
ri = (hr:hr:(Rs-hr))'; 
ri_minus_half = -1/2*hr + ri; ri_plus_half = 1/2*hr + ri;
Ri_minus_square = (ri_minus_half./ri).^2; Ri_plus_square = (ri_plus_half./ri).^2;

a=diag([-3; 1 + sigma*(Ri_minus_square+Ri_plus_square); -3])...
    +diag([4; -sigma*Ri_plus_square], 1)...
    +diag([-sigma*Ri_minus_square; 4], -1);
a(1,3) = -1; a(mr+2,mr) = -1;

b = zeros(mr+2,1); S = zeros(mr+2,1);

params.speciesAtSolid_Anode_a = a;
params.speciesAtSolid_Anode_b = b;
params.speciesAtSolid_Anode_S = S;

end