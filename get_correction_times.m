function [correction_times, delta_t] = get_correction_times(tDhat, sats, Eph)
N_sats = numel(sats);
correction_times = zeros(N_sats, 1);
delta_t = zeros(N_sats, 1);
for i=1:N_sats
    eph_col_inx = find_eph( Eph, sats(i), tDhat(i));
    eph_col = Eph(:, eph_col_inx);
    [~, tcorr] = tx_RAW2tx_GPS(tDhat(i), eph_col);
    correction_times(i) = tcorr;
    delta_t(i) = Eph(19, eph_col_inx);
end
end
