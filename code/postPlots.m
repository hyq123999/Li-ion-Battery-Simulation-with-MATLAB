%%
sz_len_cs = 7; sz_len_as = 6; sz_rates = 6; 

SOC_c_vs_len_c_rate = zeros(sz_len_cs, sz_rates); 
for i = 1:sz_len_cs
    for k = 1:sz_rates
        SOC_c_vs_len_c_rate(i, k) = SOC_c_s(i,3,k);
    end
end

SOC_a_vs_len_a_rate = zeros(sz_len_as, sz_rates); 
for j = 1:sz_len_as
    for k = 1:sz_rates
        SOC_a_vs_len_a_rate(j, k) = SOC_a_s(3,j,k);
    end
end

%%
h1 = figure;
A1 = axes;
mesh(disc_rates, len_c_ratios, SOC_c_vs_len_c_rate);
set(A1,'XScale','log');
ylabel('relative positive electrode (LiCoO_2) length'); 
xlabel('discharge rates');
zlabel('remained SOC (LiCoO_2)');
grid on
box on
view(90,90);

h2 = figure;
A2 = axes;
mesh(disc_rates, len_a_ratios, SOC_a_vs_len_a_rate);
set(A2,'XScale','log');
ylabel('relative negative electrode (LiC_6) length'); 
xlabel('discharge rates');
zlabel('remained SOC (LiC_6)');
grid on
box on
view(90,90);