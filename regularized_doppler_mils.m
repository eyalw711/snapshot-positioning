function [ellHat, bHat, resid, ns, iter_ell, rt] = regularized_doppler_mils(ellBar, presumed_time, code_phase_obs, doppler_obs, sats, Eph)
%REGULARIZED_DOPPLER_MILS This function calculates receiver position according to
% the Regularized Mixed-Integer Least-Squares approach.
%   ellBar is initial guess for receiver position (assistance)
%   presumed_time is initial guess for current GPS_Time at the receiver
%   code_phase_obs are the code phase observations relating to the sats
%   argument.
%   doppler_obs are the m/s doppler observations relating to the sats
%   argument.
%   sats is the satellite SV IDs corresponding to the measurements
%   Eph is assistance ephemeris

%%% TODO: TRY TO DETECT WHEN INTEGERS STOP UPDATING AND THEN TRY TO SOLVE
%%% BETTER WITHOUT THE REGULARIZATION.
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
bBar = 0;                       % estimated receiver-satellites clock bias
fdBar = 0;                      % frequency offset
niter = 3;                      % number of iterations
d_coarse = 76.5e-3;             % average value for signal time-of-flight satellite to earth surface
sigmaCode = 10e-9;              % code phase std
sigmaA     = 100e3/c;
sigmaD = 1800e3/c;

N_sats = numel(sats);

% Weight matrix for DOP-MILS
w = zeros(2*N_sats, 1);
w(1:N_sats,1)   = (1/sigmaCode) * ones(N_sats,1);
w(N_sats+1:2*N_sats) = (1/sigmaD   ) * ones(N_sats,1);
W = diag(w);

%%%%%%%%%%%%%
% algorithm %
%%%%%%%%%%%%%

iter_ell = ellBar;

presumed_arrival_times = presumed_time + code_phase_obs*tcode;
tDhat = presumed_arrival_times - d_coarse - bBar;   % tDhat is presumed departure times

[distances, J0, ~] = model(ellBar, tDhat, sats, Eph);
[~, range_rates, range_accs, sat_vels, es] = model_ext(ellBar, tDhat, sats, Eph);

[correction_times, ~] = get_correction_times(tDhat, sats, Eph);

% ILS method for ns calculation
A_top = [J0/c + [zeros(N_sats,d) ones(N_sats,1) ] zeros(N_sats, 1)];
B_top = -eye(N_sats)*tcode;

A_bot = zeros(N_sats, 5);
for i = 1:N_sats
    A_bot(i, 1:3) = -(1/distances(i))*(sat_vels(:,i)-es(:,i)*range_rates(i))';
    A_bot(i, 4) = range_accs(i);
    A_bot(i, 5) = range_rates(i) + c;
end
B_bot = zeros(N_sats);

rhs_top = code_phase_obs*tcode - distances/c + correction_times - bBar;
rhs_bot = -doppler_obs - range_rates - c*fdBar;

if ASSUME_HEIGHT % need at least 5 rows in the matrix for 5 variables, but can do it with 4 satellites
    A_top = [A_top; [ellBar'/norm(ellBar) 0 0]];
    B_top = [B_top; zeros(1, N_sats)];
    rhs_top = [rhs_top; 0];
    fix_W = diag(W);
    fix_W = [fix_W(1:N_sats); 1/sigmaCode; fix_W(N_sats+1:end)];
    W = diag(fix_W);
end

[x, ns] = smils(W*[A_top; A_bot] , W*[B_top; B_bot] ,W*[rhs_top; rhs_bot]);

% according to x and ns set the initial guess
ellHat = ellBar + x(1:3);
bHat = bBar + x(4);
fdHat = fdBar + x(5);

iter_ell = [iter_ell ellBar];

% ns_corr = 1;
% while sum(abs(ns_corr))
for it = 1:niter
    % Weight matrix for DOP-MILS
    w = zeros(2*N_sats, 1);
    w(1:N_sats,1)   = (1/sigmaCode) * ones(N_sats,1);
    w(N_sats+1:2*N_sats) = (1/sigmaD   ) * ones(N_sats,1);
    W = diag(w);
    
    [distances, ~, satspos] = model(ellHat, tDhat, sats, Eph); % for improving transmit times
    
    [elevs, ~] = satellite_elevations(ellHat, satspos);
    trops = tropospheric_model(elevs);
    
    tDhat = presumed_arrival_times - (ns + code_phase_obs)*tcode; %- distances/c + correction_times - bHat - trops/c;   % tDhat is presumed departure times
    correction_times = get_correction_times(tDhat, sats, Eph);
    
    [distances, J1, satspos] = model(ellHat, tDhat, sats, Eph);
    [~, range_rates, range_accs, sat_vels, es] = model_ext(ellHat, tDhat, sats, Eph);
    
    A_top = [J1/c + [zeros(N_sats,d) ones(N_sats,1) ] zeros(N_sats, 1)];
    B_top = -eye(N_sats)*tcode;
    
    satspos_rot = zeros(size(satspos));                              % earth rotation
    for s = 1:N_sats
        satspos_rot(:,s) = e_r_corr(distances(s)/c, satspos(:,s));
    end
    distances_rot = vecnorm(satspos_rot - ellHat)';                  % distance from rotated positions of satellites
    es_rot = (satspos_rot - ellHat) ./ vecnorm(satspos_rot - ellHat);
    
    A_bot = zeros(N_sats, 5);
    for i = 1:N_sats
        A_bot(i, 1:3) = -(1/distances_rot(i))*(sat_vels(:,i)-es_rot(:,i)*range_rates(i))';
        A_bot(i, 4) = range_accs(i);
        A_bot(i, 5) = range_rates(i) + c;
    end
    B_bot = zeros(N_sats);
    
    rhs_top = (ns + code_phase_obs)*tcode - distances_rot/c - trops/c + correction_times - bHat;
    rhs_bot = -doppler_obs - range_rates - fdHat;
    
    if ASSUME_HEIGHT % need at least 5 rows in the matrix for 5 variables, but can do it with 4 satellites
        A_top = [A_top; [ellBar'/norm(ellBar) 0 0]];
        B_top = [B_top; zeros(1, N_sats)];
        rhs_top = [rhs_top; 0];
        fix_W = diag(W);
        fix_W = [fix_W(1:N_sats); 1/sigmaCode; fix_W(N_sats+1:end)];
        W = diag(fix_W);
    end
    
    [x, ns_corr] = smils(W*[A_top;A_bot] , W*[B_top; B_bot] ,W*[rhs_top; rhs_bot]);

    ellHat = ellHat + x(1:3);
    bHat = bHat + x(4);
    fdHat = fdHat + x(5);
    ns = ns + ns_corr;
    
    iter_ell = [iter_ell ellHat];
end
resid = norm([rhs_top; rhs_bot]);
rt = toc;

%[ellHat, bHat, resid, ns, iter_ell2, rt_nested] = regularized_mils(ellHat, presumed_time - bHat, code_phase_obs, sats, Eph, ns);
%iter_ell = [iter_ell iter_ell2];

%rt = rt + rt_nested;
end

