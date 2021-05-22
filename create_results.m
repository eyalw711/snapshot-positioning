% create all results

% 1. CDF figure for fix accuracy with 20km 20sec error:
demo_algos;

% 2. convergence maps
basepath = pwd;
van_digglen_input = fullfile(basepath, 'results', 'final_convs', 'vd_2021_04_03__13_12_28.csv');
fernandez_input = fullfile(basepath, 'results', 'final_convs', 'shadow_d_100steps_2021_05_07__10_33_02_2021_05_07__11_09_29.csv');
mils_a_priori_input = fullfile(basepath, 'results', 'final_convs', 'apmils_2021_04_03__12_18_39.csv');
mils_dop_input = fullfile(basepath, 'results', 'final_convs', 'one_stage_reg_dop_mils_1800km_2021_05_15__13_14_57.csv');

convergence_processing({van_digglen_input, fernandez_input, mils_a_priori_input, mils_dop_input});

% 3. num satellites effect
codephases_input = fullfile(basepath, 'results', 'nsats_analysis', 'codephases_13_100steps_nsats_2021_05_15__17_14_37.csv');
dopplers_input = fullfile(basepath, 'results', 'nsats_analysis', 'dopplers_13_100steps_nsats_2021_05_15__14_10_51.csv');

convergence_processing_nsats({codephases_input, dopplers_input});
