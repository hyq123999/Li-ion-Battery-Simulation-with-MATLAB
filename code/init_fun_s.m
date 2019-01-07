function [ SOC_c, SOC_a, outputV, er ] = init_fun_s(len_c_ratio, len_a_ratio, rate)
% % clear all;
% close all;

params = getParameters();

%% user's selection
params.len_c = 80e-6 * len_c_ratio; % [m];
params.len_a = 90e-6 * len_a_ratio; % [m];
rate         = + rate; % '+': charge; '-': discharge

params.tran_num_plus = 0.364;
% params.tran_num_plus = 0.95;
% params.De = params.De / 1e2;

h = 5e-6; % spatial step [m]
hr = 5e-7; % another spatial step for solid material spheres [m]
k = 1; % time step [s]
%% more user's selection
timeLimit = 1e5; % total time [s]

ce_init = 1000; % [mol/m^3]

if (rate > 0)
    SOC_init = 0.10;
else
    SOC_init = 0.90;
end
cs_c_init = SOC_init     * params.cs_c_max;
cs_a_init = (1-SOC_init) * params.cs_a_max;

cutover_voltage = 4.3; % [V]
cutoff_voltage  = 2.5; % [V]

cutover_SOC     = 0.98;
cutoff_SOC      = 0.02; 

%% add more parameters
params = addMoreParameters(params, rate, h, hr, k);

%% create field for res
res = createEmptyResults_mat(timeLimit, k, params); %res = createEmptyResults( );

%% loop (including init)
t = 0; % start from 0
idx = 1; % start from 1
tic;
while (true)
    % j
    j_c = + params.IoverA / (params.F * params.as_c * params.len_c); % scalar
    j_a = - params.IoverA / (params.F * params.as_a * params.len_a); % scalar
    res.j_c(:,idx) = ones(params.m_c+2, 1)*j_c; %res.j_c = [res.j_c ones(params.m_c+2, 1)*j_c];
    res.j_a(:,idx) = ones(params.m_a+2, 1)*j_a; %res.j_a = [res.j_a ones(params.m_a+2, 1)*j_a];
    
    % ce
    if t == 0
        res.ce(:,idx) = ones(params.m+2, 1) * ce_init; %res.ce = ones(params.m+2, 1) * ce_init;
    else
        j = [res.j_c(:,idx);...
             zeros(params.m_s,1);...
             res.j_a(:,idx)];
        [ce, flag] = speciesAtElectrolyte(params, res.ce(:,idx-1), j); % note (idx-1), not idx
        
        if (flag == -1)
            res.exitReason = 'Electrolyte has been depleted somewhere';
            break; % exit early
        end

        res.ce(:,idx) = ce; %res.ce = [res.ce ce];

    end
    
    % cs, css
    if t == 0
        cs_c = ones(params.m_c+2, params.mr_c+2) * cs_c_init; % (m_c+2)*(mr_c+2) matrix 
        cs_a = ones(params.m_a+2, params.mr_a+2) * cs_a_init; % (m_a+2)*(mr_a+2) matrix

        res.css_c(:,idx) = ones(params.m_c+2, 1) * cs_c_init; %res.css_c = ones(params.m_c+2, 1) * cs_c_init;
        res.css_a(:,idx) = ones(params.m_a+2, 1) * cs_a_init; %res.css_a = ones(params.m_a+2, 1) * cs_a_init;
    else
        css_c = zeros(params.m_c+2, 1);
        css_a = zeros(params.m_a+2, 1);
        % % flag = 0; % important for BOTH for-loop and while-loop
        
        for i = 1:(params.m_c+2)
            if i == 1
                [cs_c_ith_new, flag] = speciesAtSolid(params, cs_c(i,:), res.j_c(i, idx), 'cathode');

                if (flag == -1)
                    res.exitReason = 'concentration of solid particles is less than 0 (somewhere at the cathode)';
                elseif (flag == 1)
                    res.exitReason = 'concentration of solid particles is greater than max value (somewhere at the cathode)';    
                end

                cs_c(i,:) = cs_c_ith_new;
                css_c(i) = cs_c_ith_new(end);
                
            else
                cs_c(i,:) = cs_c(i-1,:);
                css_c(i) = css_c(i-1);
            end
        end
        
        for i = 1:(params.m_a+2)
            if i == 1
                [cs_a_ith_new, flag] = speciesAtSolid(params, cs_a(i,:), res.j_a(i, idx), 'anode');

                if (flag == -1)
                    res.exitReason = 'concentration of solid particles is less than 0 (somewhere at the anode)';
                elseif (flag == 1)
                    res.exitReason = 'concentration of solid particles is greater than max value (somewhere at the anode)';
                end

                cs_a(i,:) = cs_a_ith_new;
                css_a(i) = cs_a_ith_new(end);
            
            else
                cs_a(i,:) = cs_a(i-1,:);
                css_a(i) = css_a(i-1);
            end               
        end
        
        res.css_c(:,idx) = css_c; %res.css_c = [res.css_c css_c];
        res.css_a(:,idx) = css_a; %res.css_a = [res.css_a css_a];
    end
    
    % SOC
    SOC_c = estimateSOC(params, cs_c, 'cathode');
    
    if (SOC_c > cutover_SOC)
        res.exitReason = 'SOC is greater than cutover SOC (at the cathode)';
    elseif (SOC_c < cutoff_SOC)
        res.exitReason = 'SOC is less than cutoff SOC (at the cathode)'; 
    end
    
    res.SOC_c(:,idx) = SOC_c; %res.SOC_c = [res.SOC_c SOC_c];
    
    SOC_a = estimateSOC(params, cs_a, 'anode');
    
    if (SOC_a > cutover_SOC)
        res.exitReason = 'SOC is greater than cutover SOC (at the anode)';
    elseif (SOC_a < cutoff_SOC)
        res.exitReason = 'SOC is less than cutoff SOC (at the anode)';
    end
    
    res.SOC_a(:,idx) = SOC_a; %res.SOC_a = [res.SOC_a SOC_a];
    
    % U
    U_c = calculateOCV(res.css_c(:,idx), params, 'cathode');
    U_a = calculateOCV(res.css_a(:,idx), params, 'anode');
    res.U_c(:,idx) = U_c; %res.U_c = [res.U_c U_c];
    res.U_a(:,idx) = U_a; %res.U_a = [res.U_a U_a];
    
    % eta
    i0_c = params.F * params.k_c * sqrt(res.ce(1:(params.m_c+2),idx)...
                                    .* (params.cs_c_max - res.css_c(:,idx))...
                                    .* res.css_c(:,idx));
    i0_a = params.F * params.k_a * sqrt(res.ce((params.m_c+params.m_s+3):end,idx)...
                                    .* (params.cs_a_max - res.css_a(:,idx))...
                                    .* res.css_a(:,idx));
    res.i0_c(:,idx) = i0_c; %res.i0_c = [res.i0_c i0_c];
    res.i0_a(:,idx) = i0_a; %res.i0_a = [res.i0_a i0_a];
    
    i0_c_bar = sum(res.i0_c(:,idx)) / numel(res.i0_c(:,idx));
    i0_a_bar = sum(res.i0_a(:,idx)) / numel(res.i0_a(:,idx));
    res.i0_c_bar(:,idx) = i0_c_bar; %res.i0_c_bar = [res.i0_c_bar i0_c_bar];
    res.i0_a_bar(:,idx) = i0_a_bar; %res.i0_a_bar = [res.i0_a_bar i0_a_bar];
    
    eta_c = (params.R*params.T/0.5/params.F) * asinh(1/2 * params.F * (+res.j_c(:,idx)) ./ res.i0_c(:,idx)); 
    eta_a = (params.R*params.T/0.5/params.F) * asinh(1/2 * params.F * (+res.j_a(:,idx)) ./ res.i0_a(:,idx)); 
    res.eta_c(:,idx) = eta_c; %res.eta_c = [res.eta_c eta_c];
    res.eta_a(:,idx) = eta_a; %res.eta_a = [res.eta_a eta_a];

    eta_c_bar = sum(res.eta_c(:,idx)) / numel(res.eta_c(:,idx)); 
    eta_a_bar = sum(res.eta_a(:,idx)) / numel(res.eta_a(:,idx));
    res.eta_c_bar(:,idx) = eta_c_bar; %res.eta_c_bar = [res.eta_c_bar eta_c_bar];
    res.eta_a_bar(:,idx) = eta_a_bar; %res.eta_a_bar = [res.eta_a_bar eta_a_bar];
    
    eta_c_bar_alt = (params.R*params.T/0.5/params.F) * asinh(1/2 * params.F * (+j_c) / res.i0_c_bar(:,idx)); 
    eta_a_bar_alt = (params.R*params.T/0.5/params.F) * asinh(1/2 * params.F * (+j_a) / res.i0_a_bar(:,idx)); 
    res.eta_c_bar_alt(:,idx) = eta_c_bar_alt; %res.eta_c_bar_alt = [res.eta_c_bar_alt eta_c_bar_alt];
    res.eta_a_bar_alt(:,idx) = eta_a_bar_alt; %res.eta_a_bar_alt = [res.eta_a_bar_alt eta_a_bar_alt];
    
    % delta_phi_e
    kappa = calculateKappa(res.ce(:,idx), params);
    res.kappa(:,idx) = kappa; %res.kappa = [res.kappa kappa];

    kappa_bar = sum(res.kappa(:,idx)) / numel(res.kappa(:,idx));
    res.kappa_bar(:,idx) = kappa_bar; %res.kappa_bar = [res.kappa_bar kappa_bar];

    phie_del_bar_sj = (params.len_c + 2*params.len_s + params.len_a) / (2*res.kappa_bar(idx))...
                            * params.IoverA;
    phie_del_bar_dm = (2*params.R*params.T/params.F) * (1-params.tran_num_plus)...
                            * 1 * (log(res.ce(1,idx)) - log(res.ce(end,idx)));
    res.phie_del_bar_sj(:,idx) = phie_del_bar_sj; %res.phie_del_bar_sj = [res.phie_del_bar_sj   phie_del_bar_sj                  ];
    res.phie_del_bar_dm(:,idx) = phie_del_bar_dm; %res.phie_del_bar_dm = [res.phie_del_bar_dm   phie_del_bar_dm                  ];
    res.phie_del_bar(:,idx)    = phie_del_bar_sj + phie_del_bar_dm; %res.phie_del_bar    = [res.phie_del_bar      phie_del_bar_sj + phie_del_bar_dm];
    
    % output voltage
    %%outputV = res.U_c(1,end) - res.U_a(end,end);
    outputV = res.U_c(1,idx) - res.U_a(end,idx)...
        + res.eta_c_bar(idx) - res.eta_a_bar(idx)...
        + res.phie_del_bar(idx);
    
    if (outputV > cutover_voltage)
        res.exitReason = 'Output voltage has reached cutover voltage';
    elseif (outputV < cutoff_voltage  && outputV ~= -Inf)
        res.exitReason = 'Output voltage has reached cutoff voltage';
    end
    
    res.outputV(:,idx) = outputV; %res.outputV = [res.outputV outputV];

    %
    if (t > timeLimit)
        res.exitReason = 'Setting time limit has been reached';
    end
    
    res.time(:,idx) = t; %res.time = [res.time t]; % log time

    % check exit and also update time
    if (strcmp(res.exitReason, ''))        
        t   = t   + k;
        idx = idx + 1;
    else
        break;  
    end
              
end
timeElapsed = toc;

%%
fprintf('Elapsed time is about %6.4f second(s).\n', timeElapsed);
fprintf('Simulated time is about %6.4f second(s).\n', res.time(idx));
fprintf(['Reason for exit: ' res.exitReason '.\n']);
fprintf(['len_c_ratio is: ' num2str(len_c_ratio) '; len_a_ratio is: ' num2str(len_a_ratio)...
         '; rate is: ' num2str(rate)  ' C.\n']);
fprintf(['SOC_c is: ' num2str(res.SOC_c(idx)) '; SOC_a is: ' num2str(res.SOC_a(idx))...
         '; output voltage is: ' num2str(res.outputV(idx))  ' V.\n']);
fprintf(' \n');

%%
SOC_c = res.SOC_c(idx); SOC_a = res.SOC_a(idx);
outputV = res.outputV(idx);

switch res.exitReason
    % ce
    case 'Electrolyte has been depleted somewhere'
        er = 11;
    % cs
    case 'concentration of solid particles is less than 0 (somewhere at the cathode)'
        er = 21;
    case 'concentration of solid particles is greater than max value (somewhere at the cathode)'
        er = 22;
    case 'concentration of solid particles is less than 0 (somewhere at the anode)'
        er = 23;
    case 'concentration of solid particles is greater than max value (somewhere at the anode)'
        er = 24;
    % soc
    case 'SOC is less than cutoff SOC (at the cathode)'
        er = 31;
    case 'SOC is greater than cutover SOC (at the cathode)'
        er = 32;
    case 'SOC is less than cutoff SOC (at the anode)'
        er = 33;
    case 'SOC is greater than cutover SOC (at the anode)'
        er = 34;
    % v  
    case 'Output voltage has reached cutoff voltage'
        er = 41;
    case 'Output voltage has reached cutover voltage'
        er = 42;
    % time  
    case 'Setting time limit has been reached'
        er = 51;
    % ??   
    otherwise
        er = -1;
end
        
        
        
        
        
        
%% trim res
res  = trimResults_mat(res, idx );

%%
% % f1 = figure;
% % set(gcf, 'Position', [100, 100, 1600, 600]);
% % time = res.time;
% % 
% % subplot(2,4,1); plot(time, res.css_c(1,:)); hold all; plot(time, res.css_a(end,:)); 
% % ylabel('Li concentration at solid particle surface {\itc}_{ss} [mol/m^3]'); xlabel('time [s]');
% % lgd = legend('cathode','anode'); lgd.Location = 'best';
% % 
% % subplot(2,4,2); plot(time, res.SOC_c); hold all; plot(time, res.SOC_a); 
% % ylabel('State of charge'); xlabel('time [s]');
% % lgd = legend('cathode','anode'); lgd.Location = 'best';
% % 
% % subplot(2,4,3); 
% % % % mesh(res.ce); view(90, 0);
% % sz = size(res.ce, 2);
% % idx = floor([1 sz/40 2*sz/40 3*sz/20 7*sz/20 15*sz/20]);
% % hold all;
% % for i = 1:numel(idx)
% %     plot(res.ce(:,idx(i)));
% % end
% % lgd = legend([num2str(res.time(idx(1))) ' s'], [num2str(res.time(idx(2))) ' s'],...
% %        [num2str(res.time(idx(3))) ' s'], [num2str(res.time(idx(4))) ' s'],...
% %        [num2str(res.time(idx(5))) ' s'], [num2str(res.time(idx(6))) ' s']);
% % lgd.Location = 'best';
% % ylabel('Li ion concentration at electrolyte {\itc}_e [mol/m^3]'); 
% % xlabel('Control volume index (cathode->separator->anode)');
% % 
% % 
% % subplot(2,4,4); 
% % % % mesh(res.kappa); view(90, 0);
% % sz = size(res.kappa, 2);
% % idx = floor([1 sz/40 2*sz/40 3*sz/20 7*sz/20 15*sz/20]);
% % hold all;
% % for i = 1:numel(idx)
% %     plot(res.kappa(:,idx(i)));
% % end
% % lgd = legend([num2str(res.time(idx(1))) ' s'], [num2str(res.time(idx(2))) ' s'],...
% %        [num2str(res.time(idx(3))) ' s'], [num2str(res.time(idx(4))) ' s'],...
% %        [num2str(res.time(idx(5))) ' s'], [num2str(res.time(idx(6))) ' s']);
% % lgd.Location = 'best';
% % ylabel('Electrolyte conductivity \kappa [S/m]');
% % xlabel('Control volume index (cathode->separator->anode)');
% % 
% % 
% % % % f2 = figure;
% % % % movegui(f2,'southwest');
% % 
% % subplot(2,4,6); 
% % plot(time, res.eta_c_bar); hold all; plot(time, res.eta_a_bar);
% % plot(time, res.eta_c_bar_alt, '--'); plot(time, res.eta_a_bar_alt, '--');
% % ylabel('Activation overpotential \eta [V]'); xlabel('time [s]');
% % lgd = legend('cathode','anode'); lgd.Location = 'best';
% % 
% % subplot(2,4,5); plot(time, res.U_c(1,:)); hold all; plot(time, res.U_a(end,:)); plot(time, res.U_c(1,:)-res.U_a(end,:)); 
% % lgd = legend('{\itU}_c','{\itU}_a','{\itU}_c-{\itU}_a'); lgd.Location = 'west';
% % ylabel('Open circuit voltage [V]'); xlabel('time [s]');
% % 
% % subplot(2,4,7);  plot(time, res.phie_del_bar); hold all; plot(time, res.phie_del_bar_sj);plot(time, res.phie_del_bar_dm);
% % lgd = legend('total','ionic flux','diff. migr.'); lgd.Location = 'east';
% % ylabel('Electrolyte potential difference along cell \Delta\phi_e [V]'); xlabel('time [s]');
% % 
% % subplot(2,4,8); plot(time, res.outputV); 
% % ylabel('Output voltage [V]'); xlabel('time [s]');
      

