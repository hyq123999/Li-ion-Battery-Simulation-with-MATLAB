function [ res ] = createEmptyResults_mat(timeLimit, k, params )

sz = ceil(timeLimit/k) + 1;

%%
    res.time = zeros(1,sz);
    res.exitReason = '';

    res.j_c = zeros(params.m_c+2,sz);
    res.j_a = zeros(params.m_a+2,sz);

    res.ce = zeros(params.m+2,sz);

    res.css_c = zeros(params.m_c+2,sz);
    res.css_a = zeros(params.m_a+2,sz);

    res.SOC_c = zeros(1,sz);
    res.SOC_a = zeros(1,sz);

    res.U_c = zeros(params.m_c+2,sz);
    res.U_a = zeros(params.m_a+2,sz);

    res.i0_c = zeros(params.m_c+2,sz);
    res.i0_a = zeros(params.m_a+2,sz);
    res.i0_c_bar = zeros(1,sz);
    res.i0_a_bar = zeros(1,sz);

    res.eta_c = zeros(params.m_c+2,sz);
    res.eta_a = zeros(params.m_a+2,sz);
    res.eta_c_bar = zeros(1,sz);
    res.eta_a_bar = zeros(1,sz);
    res.eta_c_bar_alt = zeros(1,sz);
    res.eta_a_bar_alt = zeros(1,sz);

    res.kappa = zeros(params.m+2,sz);
    res.kappa_bar = zeros(1,sz);

    res.phie_del_bar_sj = zeros(1,sz); % ionic flux as source term
    res.phie_del_bar_dm = zeros(1,sz); % diffusional migration?
    res.phie_del_bar = zeros(1,sz);

    res.outputV = zeros(1,sz);

end