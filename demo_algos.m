% This scripts generates a graph of performace evaluation of the three
% algorithms:
% - Van Diggelen's Original Shadowing Algorithm
% - MILS with A Priori Regularization
% - MILS with Doppler Regularization

clear all;
clc;

d=3;
tcode = 1e-3;

colors = [...
    0.9290, 0.6940, 0.1250;...  % Y
    0.8500, 0.3250, 0.0980;...  % R
    0, 0.4470, 0.7410;...       % B
    ];

% errors
posAssistErrorMags = [20e3 150e3];
clockBiasMags = [20 150];
assert(numel(posAssistErrorMags) == numel(clockBiasMags),...
    'must supply same error vectors of same length');

% inputs database
scenario_n_file = 'inputs/bshm1460.20n';           % may 25 2020 ephemeris
scenario_o_file = 'inputs/MEASX-VESPER-1460.txt';  % may 25 2020 ublocks observations on university roof

% ground truth
gt_ecef =  [4440082.76; 3086578.96; 3371015.61];   % receiver ground truth position
gt_lla = ecef2lla(gt_ecef');                       % receiver ground truth position in LLA format

figure; clf;

% Read RINEX ephemerides file and convert to internal Matlab format
rinexe(scenario_n_file, 'eph.dat');
Eph = get_eph('eph.dat');

% open ubx observation file
epochs = ubx_reader_codephase(scenario_o_file);
n_epochs = numel(epochs);

algo_pos_errs = zeros(n_epochs, 3);

drop_sats_vec = [-1 5 4];

subplot_inx = 0;
for err_inx = 1:numel(posAssistErrorMags)
    posAssistErrorMag = posAssistErrorMags(err_inx);
    clockBiasMag = clockBiasMags(err_inx);
    
    for drop_inx = 1:numel(drop_sats_vec)
        drop_sats = drop_sats_vec(drop_inx);
        if drop_sats < 0
            drop_sats_str = 'All Satellites In View';
        else
            drop_sats_str = num2str(drop_sats);
        end
        
        subplot_inx = subplot_inx + 1;
        subplot(numel(posAssistErrorMags), numel(drop_sats_vec), subplot_inx);
        hold on;
        set(gca, 'XScale', 'log');
        
        legend_args = {};
        
        methods = {'shadow', 'reg_mils', 'reg_dop_mils'};
        for method_inx = 1:numel(methods)
            if strcmp(methods{method_inx}, 'shadow')
                method_str = 'Van Diggelen''s';
            elseif strcmp(methods{method_inx}, 'reg_mils')
                method_str = 'MILS (A Priori Reg)';
            elseif strcmp(methods{method_inx}, 'reg_dop_mils')
                method_str = 'MILS (Doppler Reg)';
            end
            
            if drop_sats > 0
                sats_str = sprintf('Num Sats = %d', drop_sats);
            else
                sats_str = 'All Sats';
            end
            
            legend_args = [legend_args {method_str}];
            
            fixes = zeros(n_epochs, 3);
            resids = zeros(n_epochs, 1);
            
            fail_epochs = false(n_epochs, 1);
            
            rng(711); % seed randomness
            
            for ne = 1:n_epochs
                az = 360*rand;
                small_time_err = rand;
                
                % epoch coarse time inputs
                gps_time = epochs(ne).recvTOW;
                sats1 = epochs(ne).obs.sv';
                gps_codephases = epochs(ne).obs.cp';
                
                % make position assistance
                small_pos_err = rand;
                arc_len_deg = km2deg(posAssistErrorMag * 1e-3 + small_pos_err);
                assistance_lla = reckon(gt_lla(1), gt_lla(2), arc_len_deg, az);
                assistance_ecef = lla2ecef([assistance_lla gt_lla(3)]);
                ellBar = assistance_ecef';
                
                % make time assistance and altered code phases
                perb_clock_bias = clockBiasMag + small_time_err;
                presumed_time = ceil((gps_time + perb_clock_bias)/tcode)*tcode;
                axis_beta = mod(-perb_clock_bias, tcode)/tcode;
                snapshot_codephases_obs = mod(gps_codephases - axis_beta, 1);
                snapshot_doppler_obs = epochs(ne).obs.doplMS';
                
                %%%%%%%%%%%%%%%%%%%%%%%%%
                %    ILS Appproach      %
                %%%%%%%%%%%%%%%%%%%%%%%%%
                if drop_sats > 0
                    d_coarse = 76.5e-3;                                             % average value for signal time-of-flight satellite to earth surface
                    presumed_arrival_times = presumed_time + snapshot_codephases_obs*tcode;
                    tDhat = presumed_arrival_times - d_coarse;                      % tDhat is presumed departure times
                    [~, ~, satspos] = model(ellBar, tDhat, sats1, Eph);             % just for improving transmit times
                    [elevs, ~] = satellite_elevations(ellBar, satspos);
                    %                 [~, keep_sats_inx] = maxk(elevs, drop_sats);
                    keep_sats_inx = randsample(numel(sats1), drop_sats);
                    snapshot_codephases_obs = snapshot_codephases_obs(keep_sats_inx);
                    snapshot_doppler_obs = snapshot_doppler_obs(keep_sats_inx);
                    sats1 = sats1(keep_sats_inx);
                    
                elseif drop_sats == -1 
                    if abs(clockBiasMag)/60 < 60 && posAssistErrorMag < 3000e3 % filter low satellites only if the error is small
                        d_coarse = 76.5e-3;                                             % average value for signal time-of-flight satellite to earth surface
                        presumed_arrival_times = presumed_time + snapshot_codephases_obs*tcode;
                        tDhat = presumed_arrival_times - d_coarse;                      % tDhat is presumed departure times
                        [~, ~, satspos] = model(ellBar, tDhat, sats1, Eph);             % just for improving transmit times
                        [elevs, ~] = satellite_elevations(ellBar, satspos);
                        keep_sats_inx = elevs > 10;
                        snapshot_codephases_obs = snapshot_codephases_obs(keep_sats_inx);
                        snapshot_doppler_obs = snapshot_doppler_obs(keep_sats_inx);
                        sats1 = sats1(keep_sats_inx);
                    end
                end
                
                try
                    switch methods{method_inx}
                        case 'reg_dop_mils'
                            [ellHat, bHat, resid, ns, iter_ell] = regularized_doppler_mils(ellBar, presumed_time, snapshot_codephases_obs, snapshot_doppler_obs, sats1, Eph);
                        case 'shadow'
                            [ellHat, bHat, betaHat, resid, nu, iter_ell] = shadowing_ls(ellBar, presumed_time, snapshot_codephases_obs, sats1, Eph);
                        case 'reg_mils'
                            [ellHat, bHat, resid, ns, iter_ell] = regularized_mils(ellBar, presumed_time, snapshot_codephases_obs, sats1, Eph, 76*ones(size(sats1)));
                    end
                catch ME
                    if strcmp(ME.message, 'Augumented matrix is rank defficient!') || strcmp(ME.identifier,  'MATLAB:svd:matrixWithNaNInf')
                        fail_epochs(ne) = true;
                        continue;
                    else
                        rethrow(ME);
                    end
                end
                
                fixes(ne, :) = ellHat;
                resids(ne) = norm(resid);
                
                if mod(ne, 50) == 0
                    fprintf('epoch %d...\n', ne);
                end
            end
            
            abserr = vecnorm(fixes' - gt_ecef)';
            abserr = abserr(~fail_epochs);
            
            cdf = linspace(0, 1, numel(abserr));
            x = sort(abserr);
            
            [~, p90_inx] = min(abs(cdf - 0.9));
            
            if err_inx == numel(posAssistErrorMags) && drop_inx == numel(drop_sats_vec)
                semilogx(x, cdf, 'Color', colors(method_inx, :), 'LineWidth', 2, 'DisplayName', method_str); % in the end have special legend entry
            else
                semilogx(x, cdf, 'Color', colors(method_inx, :), 'LineWidth', 2);
            end
            
            if strcmp(methods{method_inx}, 'reg_dop_mils')
                h = xline(x(p90_inx), '--', 'Color', colors(method_inx, :), 'Label', sprintf('%.0f (m)', x(p90_inx))); % 'DisplayName', sprintf('%s p90', method_str),
                % the following line skip the name of the previous plot from the legend
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
            
            fprintf('num good fixes = %d out of %d\n', numel(abserr), n_epochs);
            fprintf('num fail ILS: %d (matrix rank deficient)\n', sum(fail_epochs));
        end
        
        title(sprintf('Initial Time Error = %.1f (s)\nInitial Position Error = %.1f (km)\nNum Satellites = %s',...
            clockBiasMag, posAssistErrorMag*1e-3, drop_sats_str));
        xlabel('Final Position Absolute Error (m)'); ylabel('CDF');
        xlim([0.5 1e5]); ylim([-0.1 1.1]);
        grid on;
        set(gca,'FontSize', 14);
    end
end

fig = gcf;
set(gcf, 'Position', [100, 100, 1600, 1150]);

fig.Position(3) = fig.Position(3) + 250;

Lgnd = legend('show');
Lgnd.Position(1) = 0.01;
Lgnd.Position(2) = 0.45;
sgtitle('Absolute Error CDF - Performance Comparison', 'FontSize', 18, 'fontweight', 'bold');

fprintf('done\n');


