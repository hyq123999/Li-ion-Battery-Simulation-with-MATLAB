function [ w, flag ] = speciesAtSolid ( params, cs_dist, j_Li, type )
% m: number of "interior" points
% k: time step
% cs_dist: (m+2)*1 column vector
% j_Li: scalar

%% equation
% d(cs)/dt = 1/(r^2)*d(Ds*r^2*d(cs)/dr)/dr
% Ds*d(cs)/dr = 0 @ r = 0 
% Ds*d(cs)/dr = -j @ r = R 

%% list of paramaters

% F = params.F; % Faraday constant: in C/mol
    
switch type 
    case 'cathode'         
        % Rs = params.Rs_c;
        cs_max = params.cs_c_max;
        Ds = params.Ds_c;
        % as = params.as_c;

    case 'anode' 
        % Rs = params.Rs_a;
        cs_max = params.cs_a_max;
        Ds = params.Ds_a;
        % as = params.as_a;
end

%% 
hr = params.hr; % spatial step 
% k = params.k; % time step

switch type 
    case 'cathode'         
        mr = params.mr_c;
        a  = params.speciesAtSolid_Cathode_a;
        b  = params.speciesAtSolid_Cathode_b;
        S  = params.speciesAtSolid_Cathode_S;
    case 'anode' 
        mr = params.mr_a;
        a  = params.speciesAtSolid_Anode_a;
        b  = params.speciesAtSolid_Anode_b;
        S  = params.speciesAtSolid_Anode_S;
end


%%
% h = Rs / (m+1); % spatial step
% sigma = Ds * k / (hr*hr); % must be less than 1/2

% ri = (hr:hr:(Rs-hr))'; 
% ri_minus_half = -1/2*hr + ri; ri_plus_half = 1/2*hr + ri;
% Ri_minus_square = (ri_minus_half./ri).^2; Ri_plus_square = (ri_plus_half./ri).^2;

%%
% a: size (m+2)*(m+2)
% a=diag([-3; 1 + sigma*(Ri_minus_square+Ri_plus_square); -3])...
%     +diag([4; -sigma*Ri_plus_square], 1)...
%     +diag([-sigma*Ri_minus_square; 4], -1);
% a(1,3) = -1; a(mr+2,mr) = -1;
    
% b: size (m+2)*1
% b = zeros(mr+2,1);
% b(2:end-1) = cs_dist(2:end-1);
b = cs_dist'; b(1) = 0; b(end) = 0;

% S: size (m+2)*1
% S = zeros(mr+2,1);
S(end) = (-2*hr) * (+1*j_Li/Ds); % what should be the sign of the external source j_Li?? 

w=a \ (b+S);

%%
flag = 0;
% -1
check = find(w >= 0);
if numel(check) ~= numel(w)
    flag = -1;
end
% +1
check = find(w <= cs_max);
if numel(check) ~= numel(w)
    flag = 1;
end

end

