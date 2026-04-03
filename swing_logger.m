clear;

homeDir = getenv('HOME');
addpath(genpath(fullfile(homeDir,'tailtools/swing_plotters')));
addpath(genpath(fullfile(homeDir,'tailtools/helpers')));

meta = struct();

%%% real flight data
% ac_data = readtable('~/LOGS/swing/20260219_tuning/rpy.csv');
% ac_data = readtable('~/LOGS/swing/20260219_tuning/yaw_issue.csv');

% ac_data = readtable('~/LOGS/swing/20260304_stab/01-maiden.csv');
% ac_data = readtable('~/LOGS/swing/20260304_stab/02.Kq2_Komega10.csv');
% ac_data = readtable('~/LOGS/swing/20260304_stab/08.Kq2.5Komega12.csv');

% ac_data = readtable('~/LOGS/swing/20260306_more_stab/01.rigid_structure.csv');
% ac_data = readtable('~/LOGS/swing/20260306_more_stab/02.coupling_issue_1st_n_last_pitch.csv');
% ac_data = readtable('~/LOGS/swing/20260306_more_stab/03.training_dataset.csv');
% ac_data = readtable('~/LOGS/swing/20260306_more_stab/05.diff_eff_values.csv');

% ac_data = readtable('~/LOGS/swing/20260309_again_stab/01.wls.csv');
% ac_data = readtable('~/LOGS/swing/20260309_again_stab/09.final_0.csv');
% ac_data = readtable('~/LOGS/swing/20260309_again_stab/10.final_0_with_WLS.csv');
% ac_data = readtable('~/LOGS/swing/20260309_again_stab/11.final_0_again.csv');

% ac_data = readtable('~/LOGS/swing/20260310_optitrack/02.gamma_zoo.csv'); meta = read_meta_swing("~/LOGS/swing/20260310_optitrack/02.gamma_zoo.csv");
% ac_data = readtable('~/LOGS/swing/20260310_optitrack/03.training.csv');

% ac_data = readtable('~/LOGS/swing/20260311_training_data/01.training.csv');
% ac_data = readtable('~/LOGS/swing/20260311_training_data/03.training.csv');
% ac_data = readtable('~/LOGS/swing/20260311_training_data/05.pitch90.csv');

% ac_data = readtable('~/LOGS/swing/20260312_vbatt/06.csv');
% ac_data = readtable('~/LOGS/swing/20260312_vbatt/07.training.csv');
% ac_data = readtable('~/LOGS/swing/20260312_vbatt/08.training_plus_vbattcode.csv');

% ac_data = readtable('~/LOGS/swing/20260317_guidance/02.checks_and_maiden.csv');
% ac_data = readtable('~/LOGS/swing/20260317_guidance/04.training_CT_3.5.csv');
% ac_data = readtable('~/LOGS/swing/20260317_guidance/05.bigger_bounds_g_half.csv');

% ac_data = readtable('~/LOGS/swing/20260318_quintic/02.maiden_quintic.csv');
% ac_data = readtable('~/LOGS/swing/20260318_quintic/03.Kv_1.csv');
% ac_data = readtable('~/LOGS/swing/20260318_quintic/04.Kv_2.csv');

% ac_data = readtable('~/LOGS/swing/20260319_quintic_again/04.guid_test.csv');
% ac_data = readtable('~/LOGS/swing/20260319_quintic_again/13.csv');

% ac_data = readtable('~/LOGS/swing/20260320_synching/01.manual-guided.csv');
% ac_data = readtable('~/LOGS/swing/20260320_synching/06.guided_again.csv');
% ac_data = readtable('~/LOGS/swing/20260320_synching/08.cx_cz_good.csv');
% ac_data = readtable('~/LOGS/swing/20260320_synching/09.went_crazy.csv');

% ac_data = readtable('~/LOGS/swing/20260323_monday_tests/01.bench_test_reference_signals.csv');

% ac_data = readtable('~/LOGS/swing/20260324_guidance/09.csv');
% ac_data = readtable('~/LOGS/swing/20260324_guidance/12.uncoord_circle.csv');
% ac_data = readtable('~/LOGS/swing/20260324_guidance/13.coord_circle.csv');
% ac_data = readtable('~/LOGS/swing/20260324_guidance/16.csv');

% ac_data = readtable('~/LOGS/swing/20260326_logger_issue/01.old.csv');
% ac_data = readtable('~/LOGS/swing/20260326_logger_issue/02.newcode_erik_test.csv');
% ac_data = readtable('~/LOGS/swing/20260326_logger_issue/03.old_code_till_full_interall_000.csv');
% ac_data = readtable('~/LOGS/swing/20260326_logger_issue/04.old_code_in_data_log_dir.csv');
% ac_data = readtable('~/LOGS/swing/20260326_logger_issue/05.new_code.csv');
% ac_data = readtable('~/LOGS/swing/20260326_logger_issue/06.old_code_less_data.csv');
% ac_data = readtable('~/LOGS/swing/20260326_logger_issue/07.old_code_full_data_100Hz.csv');
% ac_data = readtable('~/LOGS/swing/20260326_logger_issue/08.datalog_100Hz.csv');

% ac_data = readtable('~/LOGS/swing/20260326_flight_tests/05.wc20_stock_mot_KK_props.csv');
% ac_data = readtable('~/LOGS/swing/20260326_flight_tests/07.maiden_75c.csv');
% ac_data = readtable('~/LOGS/swing/20260326_flight_tests/08.plus_optitrack.csv');
% ac_data = readtable('~/LOGS/swing/20260326_flight_tests/09.after_training.csv');
% ac_data = readtable('~/LOGS/swing/20260326_flight_tests/10.circle.csv');

% ac_data = readtable('~/LOGS/swing/20260327_new_mot_prop_setup/01.wc19.csv');
% ac_data = readtable('~/LOGS/swing/20260327_new_mot_prop_setup/02.wc14.csv');
% ac_data = readtable('~/LOGS/swing/20260327_new_mot_prop_setup/03.wc11.csv');
% ac_data = readtable('~/LOGS/swing/20260327_new_mot_prop_setup/08.Kq2.5_Komega10.csv');
% ac_data = readtable('~/LOGS/swing/20260327_new_mot_prop_setup/09.circle_vs3.csv');
% ac_data = readtable('~/LOGS/swing/20260327_new_mot_prop_setup/10.cx_0.3.csv');
ac_data = readtable('~/LOGS/swing/20260327_new_mot_prop_setup/11.cx_0.172.csv');
% ac_data = readtable('~/LOGS/swing/20260327_new_mot_prop_setup/12.kp1.5_kv4.csv');
% ac_data = readtable('~/LOGS/swing/20260327_new_mot_prop_setup/15.kp1.5_kv4_vs4.csv');
% ac_data = readtable('~/LOGS/swing/20260327_new_mot_prop_setup/17.vs5.csv');
% ac_data = readtable('~/LOGS/swing/20260327_new_mot_prop_setup/19.vs5.csv');

% ac_data = readtable('~/LOGS/swing/20260331_coord/01.bench_test_log_format.csv');
% ac_data = readtable('~/LOGS/swing/20260331_coord/03.check_coord_signals_again.csv');
% ac_data = readtable('~/LOGS/swing/20260331_coord/04.full_coord-uncoord_mixing.csv');
% ac_data = readtable('~/LOGS/swing/20260331_coord/05.no_wref.csv');
% ac_data = readtable('~/LOGS/swing/20260331_coord/06.Kp1_Kv3.csv');
% ac_data = readtable('~/LOGS/swing/20260331_coord/07.Kp1_kv2.csv');
% ac_data = readtable('~/LOGS/swing/20260331_coord/11.loop_vs2.csv');
% ac_data = readtable('~/LOGS/swing/20260331_coord/12.manual_ceiling.csv');

% ac_data = readtable('~/LOGS/swing/20260401_immelmann/01.test_state_estimate.csv');
% ac_data = readtable('~/LOGS/swing/20260401_immelmann/07.immelmann_longer.csv');
% ac_data = readtable('~/LOGS/swing/20260401_immelmann/08.bug_fixed.csv');
% ac_data = readtable('~/LOGS/swing/20260401_immelmann/09.immelmann_vs3.5.csv');
% ac_data = readtable('~/LOGS/swing/20260401_immelmann/10.vs3.csv');
% ac_data = readtable('~/LOGS/swing/20260401_immelmann/11.again_vs3.csv');

% ac_data = readtable('~/LOGS/swing/20260402_bugs/02.again_vs3_circle_new_eff.csv'); meta = read_meta_swing('~/LOGS/swing/20260402_bugs/02.again_vs3_circle_new_eff.csv');
% ac_data = readtable('~/LOGS/swing/20260402_bugs/03.circle_vs4.csv');
% ac_data = readtable('~/LOGS/swing/20260402_bugs/04.circle_vs5.csv');
% ac_data = readtable('~/LOGS/swing/20260402_bugs/06.back_movement_and_vs2_immelmann.csv');

%%% nps
% ac_data = readtable('~/LOGS/swing/nps/test_51.csv', 'CommentStyle', '#');

% remove duplicates
[~, keep_idx] = unique(ac_data.timestamp, 'stable');
ac_data = ac_data(keep_idx, :);

% remove last line
ac_data(end, :) = [];
ac_data.timestamp = ac_data.timestamp - ac_data.timestamp(1);

Ts = median(diff(ac_data.timestamp)); 

disp("meta:"); disp(meta);
disp("Ts:"); disp(Ts);

%% Visualization
vis_swing(ac_data, 50, true);

%%
figure('Name', 'Modes');
plot(ac_data.timestamp, ac_data.guided);

%%
figure('Name', 'Translational');
plot_translational(ac_data);

%%
figure('Name', 'Reference trajectory');
plot_ref_traj(ac_data);

%%
figure('Name', 'Body Translational');
plot_body_translational(ac_data);

%%
figure('Name', 'f cmd');
plot_spec_force_cmd(ac_data);

%%
figure('Name', 'f model');
C_X = -0.900; C_Z = -0.170; C_T_v = -0.460e-8;
plot_f_model(ac_data, C_X, C_Z, C_T_v);

%%
figure('Name', '|f_τ| SP');
plot(ac_data.timestamp, -ac_data.spec_thrust_sp);
xlabel('time [s]'); ylabel('|f_τ| SP [m/s^2]'); ylim([0, 2*9.81]);

%%
figure('Name','Sign test');
plot(ac_data.timestamp, ac_data.sign_test);

%%
figure('Name', 'Eulers ZXY order');
plot_eulers(ac_data);

%%
figure('Name', 'Rates');
plot_rates(ac_data);

%%
figure('Name', 'Angular Acceleration');
plot_ang_accel(ac_data);

%%
figure('Name', 'Actuators');
plot_actuators(ac_data);

%%
figure('Name', 'Timings');
plot_timings(ac_data);

%%
figure('Name','Voltage');
plot(ac_data.timestamp, ac_data.voltage);

%%
figure('Name','IMU accel'); 
hold on;
plot(ac_data.timestamp, ac_data.acc_x);
plot(ac_data.timestamp, ac_data.acc_y);
plot(ac_data.timestamp, ac_data.acc_z);

%% for old INDI
% figure('Name', 'Stabilization Loop');
% plot_stab_attitude(ac_data);