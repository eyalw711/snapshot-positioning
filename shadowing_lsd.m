function [ellHat, bHat, resid, iter_ell, rt] = shadowing_lsd(ellBar, presumed_time, code_phase_obs, doppler_obs, sats, Eph)
%SHADOWING_LSD Navigation based on the equation:
% deltaD = [ e_i_dot 1 -rho_i_dot_dit]*[delta_x delta_y delta_z delta_fd delta_tc]^T + eps
%   ellBar is initial guess for receiver position (assistance)
%   presumed_time is initial guess for current GPS_Time at the receiver
%   code_phase_obs are the code phase observations relating to the sats
%   argument.
%   sats is the satellite SV IDs corresponding to the measurements
%   Eph is assistance ephemeris

tic

N_UNKNOWNS = 5;
FILTER_LOW_SATS = 1; % TODO!!!

% if less than 4 satellites no hope real hope with height pseudo-meas
if numel(sats) < N_UNKNOWNS - 1
    error('at least %d observations required for shadowing_iterated_least_squares approach', N_UNKNOWNS);
end

% to compare 4 sats scenarios, otherwise compare algorithms without height
% pesudo-measurements.
ASSUME_HEIGHT = (numel(sats) < 5);

%%%%%%%%%%%%%
% constants %
%%%%%%%%%%%%%

tcode = 1e-3;                   % code length in seconds
c = physconst('LightSpeed');    % speed of light
d=3;                            % dimensionality
tcBar = 0;                       % estimated receiver-satellites clock bias
fdBar = 0;                      % frequency offset
niter = 3;                      % number of iterations
d_coarse = 76.5e-3;             % average value for signal time-of-flight satellite to earth surface
sigmaCode = 10e-9;              % code phase std
sigmaA     = 100e3/c;
sigmaD = 20e3/c;

N_sats = numel(sats);

%%%%%%%%%%%%%
% algorithm %
%%%%%%%%%%%%%

iter_ell = ellBar;

presumed_arrival_times = presumed_time + code_phase_obs*tcode;
tDhat = presumed_arrival_times - d_coarse - tcBar;   % tDhat is presumed departure times

[distances, ~, ~] = model(ellBar, tDhat, sats, Eph);
[~, range_rates, range_accs, sat_vels, es] = model_ext(ellBar, tDhat, sats, Eph);

[correction_times, ~] = get_correction_times(tDhat, sats, Eph);

% FH method

H = zeros(N_sats, 5);
for i = 1:N_sats
    H(i, 1:3) = -(1/distances(i))*(sat_vels(:,i)-es(:,i)*range_rates(i))'; % e_i_dot
    H(i, 4) = 1;
    H(i, 5) = -range_accs(i);
end


% doppler predictions:
% goes by D_i = -1/c d/dt||ell - rho_i||f_0 + fd
% so actually d/dt||ell-rho_i|| = -c/f_0 * (D_i - fd)

observed = -doppler_obs;
computed = range_rates + fdBar;

x = H\(observed - computed);

% according to x set the initial guess
ellHat = ellBar + x(1:3);
fdHat = fdBar + x(4);
tcHat = tcBar + x(5);

iter_ell = [iter_ell ellBar];

for it = 1:niter
    
    [distances, ~, satspos] = model(ellHat, tDhat, sats, Eph); % for improving transmit times
    
    [elevs, ~] = satellite_elevations(ellHat, satspos);
    trops = tropospheric_model(elevs);
    
    tDhat = presumed_arrival_times - code_phase_obs*tcode + tcHat;  % tDhat is presumed departure times
    correction_times = get_correction_times(tDhat, sats, Eph);
    
    [distances, J1, satspos] = model(ellHat, tDhat, sats, Eph);
    [~, range_rates, range_accs, sat_vels, es] = model_ext(ellHat, tDhat, sats, Eph);
    
    H = zeros(N_sats, 5);
    for i = 1:N_sats
        H(i, 1:3) = -(1/distances(i))*(sat_vels(:,i)-es(:,i)*range_rates(i))'; % e_i_dot
        H(i, 4) = 1;
        H(i, 5) = -range_accs(i);
    end
    
    observed = -doppler_obs;
    computed = range_rates + fdBar;

    x = H\(observed - computed);

    % according to x set the initial guess
    ellHat = ellHat + x(1:3);
    fdHat = fdHat + x(4);
    tcHat = tcHat + x(5);
    
    iter_ell = [iter_ell ellHat];
end
resid = norm(observed - computed);
rt = toc;
bHat = tcHat;

end

