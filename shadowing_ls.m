function [ellHat, bHat, betaHat, resid, nu, iter_ell, rt] = shadowing_ls(ellBar, presumed_time, code_phase_obs, sats, Eph)
%SHADOWING_LS This function calculates receiver
%position according to the Van Diggelen's algorithm.
%   ellBar is initial guess for receiver position (assistance)
%   presumed_time is initial guess for current GPS_Time at the receiver
%   code_phase_obs are the code phase observations relating to the sats
%   argument.
%   sats is the satellite SV IDs corresponding to the measurements
%   Eph is assistance ephemeris

rt = struct();
rt.total_elapsed = 0;
rt.count = 0;

N_UNKNOWNS = 5;
FILTER_LOW_SATS = 1;

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
bBar = 0;                       % estimated receiver-satellites clock bias
niter = 24;                     % number of iterations
d_coarse = 76.5e-3;             % average value for signal time-of-flight satellite to earth surface
sigmaCode = 10e-9;              % code phase std

%%%%%%%%%%%%%
% algorithm %
%%%%%%%%%%%%%
iter_ell = [ellBar];

N_sats = numel(sats);
presumed_arrival_times = presumed_time + code_phase_obs*tcode;
tDhat = presumed_arrival_times - d_coarse - bBar;   % tDhat is presumed departure times

[distances, J0] = model(ellBar, tDhat, sats, Eph);
[correction_times, ~] = get_correction_times(tDhat, sats, Eph);

%%% step 1 - nu determination
[~,j] = min(distances); %choose closest satellite which also should be highest

nu = zeros(N_sats, 1);
nu(j) = ceil( (distances(j)/c - code_phase_obs(j)*tcode - correction_times(j)) / tcode);
beta  = (nu(j)+code_phase_obs(j))*tcode - distances(j)/c + correction_times(j);
nu    = round( (distances/c - code_phase_obs*tcode + beta - correction_times) / tcode );

%%% step 2 - optmizing location
resid = ((nu+code_phase_obs)*tcode*c - distances + correction_times*c);
H = [ J0 ones(N_sats,1) ];

if ASSUME_HEIGHT
    height_pseudo_meas(ellBar);
end

% residMag(1) = norm(resid);
tic;

delta = H \ resid;

rt.total_elapsed = rt.total_elapsed + toc;
rt.count = rt.count + 1;

ellHat = ellBar + delta(1:d,1);
bHat = bBar + delta(d+1);  % [sec] shadowed - affects transmit time
betaHat = delta(d+2);      % [m] additive (outside) in meters

iter_ell = [iter_ell ellHat];

% improve further
for it = 1:niter
    [distances, ~, satspos] = model(ellHat, tDhat, sats, Eph); % for improving transmit times
    
    [elevs, ~] = satellite_elevations(ellHat, satspos);
    low_sats = elevs < 10;
    trops = tropospheric_model(elevs);
    
    tDhat = presumed_arrival_times - distances/c - bHat - trops/c;   % improved transmit times
    correction_times = get_correction_times(tDhat, sats, Eph);       % improved correction times
    
    [distances, J1, satspos] = model(ellHat, tDhat, sats, Eph);      % improved Jacobian
    
    satspos_rot = zeros(size(satspos));                              % earth rotation
    for s = 1:N_sats
        satspos_rot(:,s) = e_r_corr(distances(s)/c, satspos(:,s));
    end
    distances_rot = vecnorm(satspos_rot - ellHat)';                  % distance from rotated positions of satellites
    
    observed = (nu+code_phase_obs)*tcode*c;
    computed = distances_rot + betaHat + trops - correction_times*c;
    resid = observed - computed;
    H = [ J1 ones(N_sats,1) ];
    
    if FILTER_LOW_SATS && (N_sats - sum(low_sats) >= N_UNKNOWNS)
        resid(low_sats) = []; % filter low satellites when solving
        H(low_sats, :) =  []; % filter low satellites when solving
    end
    
    if ASSUME_HEIGHT
        height_pseudo_meas(ellHat);
    end
    
    tic;
    
    delta = H \ resid;
    
    rt.total_elapsed = rt.total_elapsed + toc;
    rt.count = rt.count + 1;
    
    ellHat = ellHat + delta(1:d,1);
    bHat = bHat + delta(d+1);         % [sec] shadowed - affects transmit time
    betaHat = betaHat + delta(d+2);   % [m] additive (outside) in meters
    
    iter_ell = [iter_ell ellHat];
end

resid = norm(resid);

    function height_pseudo_meas(ell)
        % this adds a pseudo-measurement and a constraint that the correction
        % is orthogonal to the normal
        positional_ellBar = ell(1:d);
        height_assumption = [positional_ellBar'/norm(positional_ellBar) 0 0];
        H = [H; height_assumption];
        resid = [resid; 0];
    end
end

