function [tx_GPS, tcorr_out] = tx_RAW2tx_GPS(tx_RAW, Eph)
% TX_GPS my function as refactoring. don't really understand yet what
% correction it does.
    t0c = Eph(21);
    dt = check_t(tx_RAW-t0c);
    tcorr = (Eph(2)*dt + Eph(20))*dt + Eph(19);
    tx_GPS = tx_RAW-tcorr;
    dt = check_t(tx_GPS-t0c);
    tcorr = (Eph(2)*dt + Eph(20))*dt + Eph(19);
    tx_GPS = tx_RAW-tcorr;
    if nargout > 1
        tcorr_out = tcorr;
    end
end