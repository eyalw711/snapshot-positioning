function[ranges, range_rates, range_accs, sat_vels, es] = model_ext(ell, td, sats, Eph)
N_sats = numel(sats);
ranges = zeros(N_sats, 1);
range_rates = zeros(N_sats, 1);
range_accs = zeros(N_sats, 1);
sat_vels = zeros(3, N_sats);
es = zeros(3, N_sats);

for sat_inx = 1:N_sats
    eph_col_inx = find_eph( Eph, sats(sat_inx), td(sat_inx));
    eph_col = Eph(:, eph_col_inx);
    
    tx_GPS = tx_RAW2tx_GPS(td(sat_inx), eph_col);
    
    rho = satpos(tx_GPS, eph_col);
    
    ranges(sat_inx) = norm(ell - rho);
    
    sat_vel = satpos(tx_GPS + 0.5, eph_col) - ...
        satpos(tx_GPS - 0.5, eph_col);
    
    sat_vels(:, sat_inx) = sat_vel;
    
    e_r = (rho - ell);
    e_r = e_r / norm(e_r);
    
    es(:, sat_inx) = e_r;
    
    range_rates(sat_inx) = e_r' * sat_vel;
    
    e_t_m_half = (satpos(tx_GPS - 0.5, eph_col) - ell);
    e_t_m_half = e_t_m_half / norm(e_t_m_half);
    
    e_t_p_half = (satpos(tx_GPS + 0.5, eph_col) - ell);
    e_t_p_half = e_t_p_half / norm(e_t_p_half);
    
    rho_t_m_1 = satpos(tx_GPS - 1, eph_col);
    rho_t_p_1 = satpos(tx_GPS + 1, eph_col);
    
    range_accs(sat_inx) = -e_t_m_half'*rho_t_m_1 - e_t_p_half'*rho_t_p_1 ...
        +(e_t_m_half + e_t_p_half)'*rho;
    
end
end