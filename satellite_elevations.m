function [elevs, azimuths] = satellite_elevations(ellbar, satspos)
N_sats = size(satspos, 2);
elevs = zeros(N_sats, 1);
azimuths = zeros(N_sats, 1);
for s = 1:N_sats
    [azimuths(s), elevs(s), ~] = topocent(ellbar, satspos(:,s)-ellbar);
end
end