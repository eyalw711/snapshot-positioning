function convergence_map_data(csv_path, algo_str, n_steps, max_pos_err_str, max_time_err_str)
%CONVERGENCE_MAP runs the convergence experiment. after done, plot with
%convergence_processing.m file

rng default
tcode = 1e-3;

tic

alg_types_lib = {'shadow', 'reg_mils', 'reg_dop_mils'};

assert(ismember(algo_str, alg_types_lib), 'invalid algo_str: %s', algo_str);

alg_types = {algo_str};

n_time_steps = str2num(n_steps);
n_pos_steps = str2num(n_steps);

one_sided_time_err = linspace(0, str2num(max_time_err_str), n_time_steps);
time_err_base = [-one_sided_time_err(end:-1:2) one_sided_time_err]; % seconds
pos_err_mag_base = linspace(0, str2num(max_pos_err_str), n_pos_steps); % meters

fprintf('running algo %s with max-time-err: %f[s], max-pos-err: %f[m]\n', algo_str, str2num(max_time_err_str), str2num(max_pos_err_str));

n_exps = 8; % random azimuth, epoch, and sub-ms error
n_buffer_rows = 1e3;
col_names = {'alg', 'time', 'ecef_x', 'ecef_y', 'ecef_z', 't_err_base',...
    't_err_small', 't_err', 'pos_err_mag_base', 'pos_err_small', 'az', 't_0', 'ecef_x0',...
    'ecef_y0', 'ecef_z0', 'ell_x', 'ell_y', 'ell_z', 'bHat', 'betaHat',...
    'resid_norm'};
n_cols = numel(col_names);
n_lines_full_file = prod([numel(alg_types), numel(time_err_base), numel(pos_err_mag_base), n_exps]);

tmstr = datestr(now, 'yyyy_mm_dd__HH_MM_SS');
[path, filename, extension] = fileparts(csv_path);
out_csv = fullfile(path, [filename '_' tmstr extension]);

f = waitbar(0, 'reading input csv');
try
    if exist(csv_path, 'file')
        tbl = readtable(csv_path, 'Delimiter', ',');
    else
        tbl = array2table(zeros(0, n_cols), 'VariableNames', col_names);
    end
    writetable(tbl, out_csv);
    buffer_tbl = array2table(zeros(0, n_cols), 'VariableNames', col_names);
    n_curr = size(tbl, 1);
    last_percent = ceil(n_curr*100/n_lines_full_file);
    
    waitbar(last_percent/100, f, 'generating epochs input');
    s_input = gen_input('ubx');
    
    waitbar(last_percent/100, f, sprintf('done %d%%', last_percent));
    
    for alg_type = alg_types
        alg_type_str = alg_type{1};
        for clockBiasMag = time_err_base
            for posAssistErrorMag = pos_err_mag_base
                n_rows = size(tbl(...
                    strcmp(alg_type_str, tbl.alg) & ...
                    clockBiasMag == tbl.t_err_base & ...
                    posAssistErrorMag == tbl.pos_err_mag_base ...
                    ,:), 1);
                n_missing_rows = n_exps - n_rows;
                
                for r_inx = 1:n_missing_rows
                    cell_row = run_algo(alg_type_str, clockBiasMag, posAssistErrorMag);
                    buffer_tbl = [buffer_tbl; cell_row];
                end
                
                if size(buffer_tbl, 1) > n_buffer_rows
                    flush_results();
                    toc
                end
            end
        end
    end
    flush_results();
    close(f);
catch ME
    fprintf(ME.message);
    close(f);
    rethrow(ME);
end


    function new_row = run_algo(alg_type_str, clockBiasMag, posAssistErrorMag)
        az = 360*rand;
        small_time_err = rand;
        ne = randi(s_input.n_epochs, 1);
        
        % epoch coarse time inputs
        gps_time = s_input.epochs(ne).recvTOW;
        sats1 = s_input.epochs(ne).obs.sv';
        gps_codephases = s_input.epochs(ne).obs.cp';
        
        % make position assistance
        small_pos_err = rand;
        arc_len_deg = km2deg(posAssistErrorMag * 1e-3 + small_pos_err);
        assistance_lla = reckon(s_input.gt_lla(1), s_input.gt_lla(2), arc_len_deg, az);
        assistance_ecef = lla2ecef([assistance_lla s_input.gt_lla(3)]);
        ellBar = assistance_ecef';
        
        % make time assistance and altered code phases
        perb_clock_bias = clockBiasMag + small_time_err;
        presumed_time = ceil((gps_time + perb_clock_bias)/tcode)*tcode;
        axis_beta = mod(-perb_clock_bias, tcode)/tcode;
        snapshot_codephases_obs = mod(gps_codephases - axis_beta, 1);
        snapshot_doppler_obs = s_input.epochs(ne).obs.doplMS';
        
        try
            betaHat = 0;
            switch alg_type_str
                case 'shadow'
                    [ellHat, bHat, betaHat, resid, ~, ~] = shadowing_ls(ellBar, presumed_time, snapshot_codephases_obs, sats1, s_input.Eph);
                case 'reg_mils'
                    [ellHat, bHat, resid, ~, ~] = regularized_mils(ellBar, presumed_time, snapshot_codephases_obs, sats1, s_input.Eph, 76*ones(size(sats1)));
                case 'reg_dop_mils'
                    [ellHat, bHat, resid, ~, ~] = regularized_doppler_mils(ellBar, presumed_time, snapshot_codephases_obs, snapshot_doppler_obs, sats1, s_input.Eph);
            end
        catch algo_me
            if strcmp(algo_me.message, 'Augumented matrix is rank defficient!')
                ellHat = zeros(size(ellBar));
                bHat = 0;
                resid = inf;
            end
        end
        
%         col_names = {'alg', 'time', 'ecef_x', 'ecef_y', 'ecef_z', 't_err_base',...
%             't_err_small', 't_err', 'pos_err_mag_base', 'pos_err_small', 'az', 't_0', 'ecef_x0',...
%             'ecef_y0', 'ecef_z0', 'ell_x', 'ell_y', 'ell_z', 'bHat', 'betaHat',...
%             'resid_norm'};

        % add row
        new_row = {alg_type_str, gps_time, s_input.gt_ecef(1), s_input.gt_ecef(2), s_input.gt_ecef(3), ...
            clockBiasMag, small_time_err, perb_clock_bias, posAssistErrorMag, small_pos_err, az, presumed_time, ...
            ellBar(1), ellBar(2), ellBar(3), ellHat(1), ellHat(2), ellHat(3), bHat, betaHat, resid};
        
    end

    function flush_results
        stc = table2struct(buffer_tbl);
        fid = fopen(out_csv, 'a');
        for j = 1:size(buffer_tbl, 1)
            fprintf(fid, '%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n',...
                stc(j).alg, stc(j).time, stc(j).ecef_x, stc(j).ecef_y, stc(j).ecef_z, ...
                stc(j).t_err_base, stc(j).t_err_small, stc(j).t_err, stc(j).pos_err_mag_base, stc(j).pos_err_small, ...
                stc(j).az, stc(j).t_0, ...
                stc(j).ecef_x0, stc(j).ecef_y0, stc(j).ecef_z0, stc(j).ell_x, stc(j).ell_y, stc(j).ell_z, ...
                stc(j).bHat, stc(j).betaHat, stc(j).resid_norm);
        end
        fclose(fid);
        
        n_curr = n_curr + size(buffer_tbl, 1);
        curr_percent = ceil(n_curr*100/n_lines_full_file);
        if curr_percent > last_percent
            last_percent = curr_percent;
            waitbar(last_percent/100, f, sprintf('done %d%%', last_percent));
        end
        buffer_tbl = array2table(zeros(0, n_cols), 'VariableNames', col_names);
    end

    function s_out = gen_input(type_str)
        if strcmp(type_str, 'ubx')
            s_out.epochs = ubx_reader_codephase('inputs/MEASX-VESPER-1460.txt');
        elseif strcmp(type_str, 'tag')
            epochs = load('inputs/snaps_epochs_1460.mat');
            s_out.epochs = epochs.epochs;
        else
            error('invalid input type');
        end
        s_out.n_epochs = numel(s_out.epochs);
        
        s_out.gt_ecef = [4440082.76 3086578.96 3371015.61];
        s_out.gt_lla = ecef2lla(s_out.gt_ecef);
        s_out.iono = [6.5193D-09  2.2352D-08 -5.9605D-08 -1.1921D-07 8.6016D+04  9.8304D+04 -6.5536D+04 -5.2429D+05];
        
        rinexe('inputs/bshm1460.20n','eph.dat');
        s_out.Eph = get_eph('eph.dat');
    end
end



