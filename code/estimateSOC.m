function [ SOC ] = estimateSOC(params, cs, type)

switch type 
    case 'cathode'         
        Rs = params.Rs_c;
        cs_max = params.cs_c_max;

    case 'anode' 
        Rs = params.Rs_a;
        cs_max = params.cs_a_max;
end

x_count = size(cs, 1);
r_count = size(cs, 2);

num = 0; den = 0;
for i = 1:x_count
    for j = 1:r_count
        num = num + cs(i, j) * j^2; % spherical shape
        den = den + j^2;
    end
end

cs_bar = num / den;
SOC = cs_bar / cs_max;
    
end