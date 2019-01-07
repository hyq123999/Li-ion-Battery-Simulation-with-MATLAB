function [ res ] = trimResults_mat(res, idx )

% sz = ceil(timeLimit/k) + 1;

    res.time = res.time(:,1:idx);
    %res.exitReason = '';

    res.j_c = res.j_c(:,1:idx);
    res.j_a = res.j_a(:,1:idx);

    res.ce = res.ce(:,1:idx);

    res.css_c = res.css_c(:,1:idx);
    res.css_a = res.css_a(:,1:idx);

    res.SOC_c = res.SOC_c(:,1:idx);
    res.SOC_a = res.SOC_a(:,1:idx);

    res.U_c = res.U_c(:,1:idx);
    res.U_a = res.U_a(:,1:idx);

    res.i0_c = res.i0_c(:,1:idx);
    res.i0_a = res.i0_a(:,1:idx);
    res.i0_c_bar = res.i0_c_bar(:,1:idx);
    res.i0_a_bar = res.i0_a_bar(:,1:idx);

    res.eta_c = res.eta_c(:,1:idx);
    res.eta_a = res.eta_a(:,1:idx);
    res.eta_c_bar = res.eta_c_bar(:,1:idx);
    res.eta_a_bar = res.eta_a_bar(:,1:idx);
    res.eta_c_bar_alt = res.eta_c_bar_alt(:,1:idx);
    res.eta_a_bar_alt = res.eta_a_bar_alt(:,1:idx);

    res.kappa = res.kappa(:,1:idx);
    res.kappa_bar = res.kappa_bar(:,1:idx);

    res.phie_del_bar_sj = res.phie_del_bar_sj(:,1:idx); % ionic flux as source term
    res.phie_del_bar_dm = res.phie_del_bar_dm(:,1:idx); % diffusional migration?
    res.phie_del_bar = res.phie_del_bar(:,1:idx);

    res.outputV = res.outputV(:,1:idx);

end