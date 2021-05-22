function [ellHat, bHat, resid, ns, iter_ell, rt] = regularized_mils(ellBar, presumed_time, code_phase_obs, sats, Eph, nBar)
%REGULARIZED_MILS This function calculates receiver position according to
% the Regularized Mixed-Integer Least-Squares approach.
%   ellBar is initial guess for receiver position (assistance)
%   presumed_time is initial guess for current GPS_Time at the receiver
%   code_phase_obs are the code phase observations relating to the sats
%   argument.
%   sats is the satellite SV IDs corresponding to the measurements
%   Eph is assistance ephemeris

%%% TODO: TRY TO DETECT WHEN INTEGERS STOP UPDATING AND THEN TRY TO SOLVE
%%% BETTER WITHOUT THE REGULARIZATION.
rt = struct();
rt.total_elapsed = 0;
rt.count = 0;

N_UNKNOWNS = 4;
FILTER_LOW_SATS = 1; % TODO!!!

% if less than 4 satellites no hope real hope with height pseudo-meas
if numel(sats) < N_UNKNOWNS
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
niter = 3;                      % number of iterations
d_coarse = 76.5e-3;             % average value for signal time-of-flight satellite to earth surface
sigmaCode = 10e-9;              % code phase std
sigmaA     = 100e3/c;

N_sats = numel(sats);

% Weight matrix for MILS
w = zeros(2*N_sats, 1);
w(1:N_sats,1)   = (1/sigmaCode) * ones(N_sats,1);
w(N_sats+1:2*N_sats) = (1/sigmaA   ) * ones(N_sats,1);
W = diag(w);

%%%%%%%%%%%%%
% algorithm %
%%%%%%%%%%%%%

iter_ell = ellBar;

presumed_arrival_times = presumed_time + code_phase_obs*tcode;
tDhat = presumed_arrival_times - d_coarse - bBar;   % tDhat is presumed departure times

[distances, J0] = model(ellBar, tDhat, sats, Eph);
[correction_times, ~] = get_correction_times(tDhat, sats, Eph);

% ILS method for ns calculation
A_s = J0/c + [zeros(N_sats,d) ones(N_sats,1) ];
B_s = -eye(N_sats)*tcode;
rhs = (nBar + code_phase_obs)*tcode - distances/c + correction_times - bBar;
A = [A_s; J0/c];
B = [B_s; zeros(N_sats)];
R = [rhs; zeros(N_sats,1)];

if ASSUME_HEIGHT
    A = [A; [ellBar'/norm(ellBar) 0]];
    B = [B; zeros(1, N_sats)];
    R = [R; 0];
    fix_W = diag(W);
    fix_W = [fix_W; 1/sigmaCode];
    W = diag(fix_W);
end

tic;

[x, ns_corr] = smils(W*A , W*B ,W*R);

rt.total_elapsed = rt.total_elapsed + toc;
rt.count = rt.count + 1;

% according to x and ns set the initial guess
ellHat = ellBar + x(1:3);
bHat = bBar + x(4);
ns = nBar + ns_corr;

iter_ell = [iter_ell ellBar];

for it = 1:niter
    % Weight matrix for MILS
    w = zeros(2*N_sats, 1);
    w(1:N_sats,1)   = (1/sigmaCode) * ones(N_sats,1);
    w(N_sats+1:2*N_sats) = (1/sigmaA   ) * ones(N_sats,1);
    W = diag(w);

    [distances, ~, satspos] = model(ellHat, tDhat, sats, Eph); % for improving transmit times
    
    [elevs, ~] = satellite_elevations(ellHat, satspos);
    trops = tropospheric_model(elevs);
    
    tDhat = presumed_arrival_times - (ns + code_phase_obs)*tcode; %distances/c + correction_times - bHat - trops/c;   % tDhat is presumed departure times
    correction_times = get_correction_times(tDhat, sats, Eph);
    
    [distances, J1, satspos] = model(ellHat, tDhat, sats, Eph);
    
    A_s = J1/c + [zeros(N_sats,d) ones(N_sats,1) ];
    B_s = -eye(N_sats)*tcode;
    
    satspos_rot = zeros(size(satspos));                              % earth rotation
    for s = 1:N_sats
        satspos_rot(:,s) = e_r_corr(distances(s)/c, satspos(:,s));
    end
    distances_rot = vecnorm(satspos_rot - ellHat)';                  % distance from rotated positions of satellites

    rhs = (ns + code_phase_obs)*tcode - distances_rot/c - trops/c + correction_times - bHat;
    
    A = [A_s;J0/c];
    B = [B_s; zeros(N_sats)];
    R = [rhs; zeros(N_sats,1)];

    if ASSUME_HEIGHT
        A = [A; [ellBar'/norm(ellBar) 0]];
        B = [B; zeros(1, N_sats)];
        R = [R; 0];
        fix_W = diag(W);
        fix_W = [fix_W; 1/sigmaCode];
        W = diag(fix_W);
    end
    
    tic;
    
    [x, ns_corr] = smils(W*A , W*B ,W*R);
    
    rt.total_elapsed = rt.total_elapsed + toc;
    rt.count = rt.count + 1;

    ellHat = ellHat + x(1:3);
    bHat = bHat + x(4);
    ns = ns + ns_corr;
    
    iter_ell = [iter_ell ellHat];
end

resid = norm(R); %([A_s; J1/c]*[ellHat; bHat] + [B_s; zeros(N_sats)]*ns - [rhs; zeros(N_sats,1)]);
end

