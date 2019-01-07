function [ U ] = calculateOCV(css, params, type)
% % cs is a vector
% % U is a vector
switch type
    case 'cathode'
       
% %         theta_p  = css ./ params.cs_c_max;
% % 
% %         % Define the OCV for the positive electrode
% %         U_c  = (-4.656+88.669*theta_p.^2 - 401.119*theta_p.^4 + 342.909*theta_p.^6 - 462.471*theta_p.^8 + 433.434*theta_p.^10);
% %         U_c  = U_c./(-1+18.933*theta_p.^2-79.532*theta_p.^4+37.311*theta_p.^6-73.083*theta_p.^8+95.96*theta_p.^10);
% % 
% %         U = U_c;
        
        SOCs = 1:-0.05:0;
        Us = [4.20 4.12 4.04 4.00...
              3.94 3.87 3.81 3.77...
              3.71 3.69 3.67 3.63...
              3.62 3.60 3.57 3.55...
              3.52 3.50 3.47 3.35...
              3.03];
         
        theta_c = css ./ params.cs_c_max;
        U = interp1(SOCs, Us, theta_c, 'spline');
                  
    
    case 'anode'
        
% %         theta_n  = css ./ params.cs_a_max;
% % 
% %         % Define the OCV for the negative electrode
% %         U_a   = 0.7222 + 0.1387*theta_n + 0.029*theta_n.^0.5 - 0.0172./theta_n + 0.0019./theta_n.^1.5 + 0.2808*exp(0.9-15*theta_n)-0.7984*exp(0.4465*theta_n - 0.4108);
% % 
% %         U = U_a;

        SOCs = 0.95:-0.05:0.05;
        Us = [      0.300 0.200 0.175...
              0.140 0.125 0.123 0.120...
              0.116 0.111 0.100 0.088...
              0.081 0.081 0.080 0.079...
              0.077 0.075 0.070 0.038
             ];
         
        theta_a = css ./ params.cs_a_max;
        U = interp1(SOCs, Us, theta_a, 'spline');


    
end

end