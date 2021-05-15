function convergence_processing_nsats(files_cell)
% this will plot the convergence maps of the paths in files_cell

colors = [...
    0.9290, 0.6940, 0.1250;...  % Y
    0.4660 0.6740 0.1880;...    % G
    0.8500, 0.3250, 0.0980;...  % R
    0, 0.4470, 0.7410;...       % B
    ];

for file_inx = 1:numel(files_cell)
    file_path = files_cell{file_inx};
    t = readtable(file_path);
    
    t.pos_err = sqrt(sum((table2array(t(:, {'ecef_x', 'ecef_y', 'ecef_z'})) - table2array(t(:, {'ell_x', 'ell_y', 'ell_z'}))).^2, 2));
    t.t_err_base = abs(t.t_err_base);
    [group, id_alg, id_nsats, id_time_err_base, id_pos_err_base] = findgroups(t.alg, t.nsats, t.t_err_base ,t.pos_err_mag_base);
    
    num_succ_exper = splitapply(@(e) sum(e<1e3) , t.pos_err, group);
    num_tot_exper = splitapply(@(e) numel(e) , t.pos_err, group);
    
    [alg_code, alg_ids] = findgroups(id_alg);
    [n_sats_code, n_sats_ids] = findgroups(id_nsats);
    [time_err_code, time_err_ids] = findgroups(id_time_err_base);
    [pos_err_code, pos_err_ids] = findgroups(id_pos_err_base);
    
    n_groups = numel(num_succ_exper);
    n_algs = numel(alg_ids);
    n_nsats = numel(n_sats_ids);
    n_time_errs = numel(time_err_ids);
    n_pos_err_ids = numel(pos_err_ids);
    
    succ_exper_map = zeros(n_algs, n_nsats, n_time_errs, n_pos_err_ids);
    tot_exper_map = zeros(n_algs, n_nsats, n_time_errs, n_pos_err_ids);
    
    for err_inx = 1:n_groups
        succ_exper_map(alg_code(err_inx), n_sats_code(err_inx), time_err_code(err_inx), pos_err_code(err_inx)) = num_succ_exper(err_inx);
        tot_exper_map(alg_code(err_inx), n_sats_code(err_inx), time_err_code(err_inx), pos_err_code(err_inx)) = num_tot_exper(err_inx);
    end
    
    portion_per_n_sats = zeros(n_nsats, n_algs);
    for alg_inx = 1:n_algs
        for sats_inx = 1:n_nsats
            curr_portion_map = squeeze(succ_exper_map(alg_inx, sats_inx, :, :));
            curr_tot_map = squeeze(tot_exper_map(alg_inx, sats_inx, :, :));
            portion_per_n_sats(sats_inx, alg_inx) = sum(curr_portion_map(:)) / sum(curr_tot_map(:));
        end
    end
    
    figure(1); hold on;
    for alg_inx = 1:n_algs
        
        if strcmp(alg_ids{alg_inx},'shadow')
            new_alg_id = 'Van Diggelen''s (Original)';
            pos_err_div = 1000;
            pos_err_str = 'km';
            tim_err_div = 1;
            tim_err_str = 's';
            clr = colors(1, :);
        elseif strcmp(alg_ids{alg_inx},'shadow_d')
            new_alg_id = 'Fernandez Hernandez Original';
            pos_err_div = 1000*1000;
            pos_err_str = 'x1000 km';
            tim_err_div = 3600;
            tim_err_str = 'Hrs';
            clr = colors(2, :);
        elseif strcmp(alg_ids{alg_inx},'reg_mils')
            new_alg_id = 'MILS with A Priori Regularization';
            pos_err_div = 1000;
            pos_err_str = 'km';
            tim_err_div = 1;
            tim_err_str = 's';
            clr = colors(3, :);
        elseif strcmp(alg_ids{alg_inx},'reg_dop_mils')
            new_alg_id = 'MILS with Doppler Regularization';
            pos_err_div = 1000*1000;
            pos_err_str = 'x1000 km';
            tim_err_div = 3600;
            tim_err_str = 'Hrs';
            clr = colors(4, :);
        end
        
        for jj = 1:n_nsats
            figure((alg_inx - 1)*n_nsats + jj + 1);
            curr_portion_map = squeeze(succ_exper_map(alg_inx, jj, :, :));
            imagesc(pos_err_ids/pos_err_div, time_err_ids/tim_err_div, curr_portion_map);
            title(sprintf('%s - %d sats', new_alg_id, n_sats_ids(jj)));
            xlabel(sprintf('initial position error (%s)', pos_err_str));
            ylabel(sprintf('initial time error (%s)', tim_err_str));
            colorbar;
            set(gca,'FontSize',14);
        end
        figure(1); hold on;
        plot(n_sats_ids, 100*portion_per_n_sats(:, alg_inx), ...
            'Color', clr, 'LineWidth', 2, 'DisplayName', new_alg_id)
    end
    
    title('Portion low error experiments VS num observations');
    xlabel('Num Observations');
    ylabel('Percent of Experiments');
    Lgnd = legend('show');
    Lgnd.Position(1) = 0.15;
    Lgnd.Position(2) = 0.850;
    
    grid on;
    set(gca,'FontSize', 14);
    set(gcf, 'Position', [100, 100, 1600, 1150]);
end
end