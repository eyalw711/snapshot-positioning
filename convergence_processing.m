function convergence_processing(files_cell)
% this will plot the convergence maps of the paths in files_cell

for file_inx = 1:numel(files_cell)
    file_path = files_cell{file_inx};
    t = readtable(file_path);
    
    t.pos_err = sqrt(sum((table2array(t(:, {'ecef_x', 'ecef_y', 'ecef_z'})) - table2array(t(:, {'ell_x', 'ell_y', 'ell_z'}))).^2, 2));
    t.t_err_base = abs(t.t_err_base);
    [group, id_alg, id_time_err_base, id_pos_err_base] = findgroups(t.alg, t.t_err_base ,t.pos_err_mag_base);
    
    res_pos_err = splitapply(@(e) sum(e<1e3)/numel(e) , t.pos_err, group);
    
    [alg_code, alg_ids] = findgroups(id_alg);
    [time_err_code, time_err_ids] = findgroups(id_time_err_base);
    [pos_err_code, pos_err_ids] = findgroups(id_pos_err_base);
    
    n_groups = numel(res_pos_err);
    n_algs = numel(alg_ids);
    n_time_errs = numel(time_err_ids);
    n_pos_err_ids = numel(pos_err_ids);
    
    conv_map = zeros(n_algs, n_time_errs, n_pos_err_ids);
    
    for err_inx = 1:n_groups
        conv_map(alg_code(err_inx), time_err_code(err_inx), pos_err_code(err_inx)) = res_pos_err(err_inx);
    end
    
    for alg_inx = 1:n_algs
        if strcmp(alg_ids{alg_inx},'shadow')
            new_alg_ids{alg_inx} = 'Van Diggelen''s (Original)';
            pos_err_div = 1000;
            pos_err_str = 'km';
            tim_err_div = 1;
            tim_err_str = 's';
        elseif strcmp(alg_ids{alg_inx},'shadow-d')
            new_alg_ids{alg_inx} = 'Fernandez Hernandez Original';
            pos_err_div = 1000*1000;
            pos_err_str = 'x1000 km';
            tim_err_div = 3600;
            tim_err_str = 'Hrs';
        elseif strcmp(alg_ids{alg_inx},'reg_mils')
            new_alg_ids{alg_inx} = 'MILS with A Priori Regularization';
            pos_err_div = 1000;
            pos_err_str = 'km';
            tim_err_div = 1;
            tim_err_str = 's';
        elseif strcmp(alg_ids{alg_inx},'reg_dop_mils')
            new_alg_ids{alg_inx} = 'MILS with Doppler Regularization';
            pos_err_div = 1000*1000;
            pos_err_str = 'x1000 km';
            tim_err_div = 3600;
            tim_err_str = 'Hrs';
        end
        
        figure;
        curr_conv_map = squeeze(conv_map(alg_inx, :, :));
        imagesc(pos_err_ids/pos_err_div, time_err_ids/tim_err_div, curr_conv_map);
        title(new_alg_ids{alg_inx});
        xlabel(sprintf('initial position error (%s)', pos_err_str));
        ylabel(sprintf('initial time error (%s)', tim_err_str));
        colorbar;
        set(gca,'FontSize',14);
    end
    set(gcf, 'Position', [100, 100, 600, 450]);
end
end