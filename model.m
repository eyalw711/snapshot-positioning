function[distances, J, satspos] = model(ell, td, sats, Eph)
N_sats = numel(sats);
J = zeros(N_sats, 4);
distances = zeros(N_sats, 1);

for r = 1:N_sats
    eph_col_inx = find_eph( Eph, sats(r), td(r));
    eph_col = Eph(:, eph_col_inx);
    
    tx_GPS = tx_RAW2tx_GPS(td(r), eph_col);
    
    rho = satpos(tx_GPS, eph_col);
    if nargout > 2
        satspos(:,r) = rho;
    end
    
    distances(r) = norm(ell - rho);
    
    sat_vel_approx = satpos(tx_GPS + 0.5, eph_col) - ...
        satpos(tx_GPS - 0.5, eph_col);
    
    e_r = (ell - rho)';
    e_r = e_r / norm(e_r);
    J(r, 1:3) = e_r;
    J(r, 4)   = e_r * sat_vel_approx;
end
end