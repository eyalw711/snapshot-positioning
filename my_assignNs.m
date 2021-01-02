function [N, N0_inx] = my_assignNs(sats, svInxListByDistance, obs, Eph, TOW_assist ,rec_loc_assist, approx_distances)
% MY_ASSIGNNS function to assign Ns according to Van-Diggelen's algorithm
    light_ms = 299792458 * 0.001;
    N = zeros(size(sats));
    approx_distances = approx_distances / light_ms; % distances in millisec
    
    N0_inx = svInxListByDistance(1);
    N(N0_inx) = floor(approx_distances(N0_inx));
    
    delta_t = Eph(19, :)' * 1000; % from sec to millisec
    
    for k = svInxListByDistance(2:end)
        N(k) = round(N(N0_inx) + obs(N0_inx) - obs(k) +...
            (approx_distances(k) - delta_t(k)) - (approx_distances(N0_inx) - delta_t(N0_inx)));
    end
end